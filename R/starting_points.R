
extract.ind = function(seq) {
     ## finds the matching indices of a sequence.
     counter = 0
     ind=NULL
     for (i in 1:length(seq)) {
          if (seq[i] == '-') {
               ind[i] = 0
          } else {
               counter = counter + 1
               ind[i] = counter
          }
     }
     return(ind)
}

SPchoose = function(data1, data2, k = 4,cut = 0.2) {
     ## Returns a matrix of starting points.
     seq1 = data1[, 1]
     seq2 = data2[, 1]
     Lmin = min(length(seq1), length(seq2))
     seq1 = paste(seq1, collapse = '')
     seq2 = paste(seq2, collapse = '')
     x0 = AAStringSet(c(seq1, seq2))
     #res1 = msa(x0,method = 'ClustalW')
     #res2 = msa(x0,method = 'ClustalOmega')
     res1 = msa(x0, method = 'Muscle', type = 'protein')
     s1 = strsplit(as.character(slot(res1, 'unmasked')[1]), '')[[1]]
     s2 = strsplit(as.character(slot(res1, 'unmasked')[2]), '')[[1]]
     ind1 = extract.ind(s1)
     ind2 = extract.ind(s2)
     ali.mat = cbind(ind1, ind2)
     gaps = which(ali.mat == 0, arr.ind = TRUE)
     ali.mat = ali.mat[-gaps[, 1], ]
     reps = 1
     d0 = 1.24 * (Lmin - 15)^(1 / 3) - 1.8
     cond = TRUE
     while (reps < 30 && cond) {
          data.proc = array(NA,dim = c(dim(ali.mat)[1], 3, 2))
          data.proc[, , 1] = as.matrix(data1[ali.mat[, 1], -1])
          data.proc[, , 2] = as.matrix(data2[ali.mat[, 2], -1])
          proc = gpa_C(data.proc)
          euc.dist = NULL
          for (i in 1:dim(ali.mat)[1]) {
               euc.dist[i] = dist(rbind(proc$rot[i, , 1],
                                        proc$rot[i, , 2]))
          }
          TM = 1 / (1 + (euc.dist^2) / (d0^2))
          #m = median(TM)
          #std = sd(TM)
          #cut = m+sig*std
          rows = which(TM < cut)
          if (length(rows) == 0) {
               break
          }
          ali.mat = ali.mat[-rows, ]
          reps = reps + 1
          if (is.null(dim(ali.mat)[1])) {
               cond = FALSE
          } else {
               if (dim(ali.mat)[1] < 3){
                    cond = FALSE
               }
          }
          
     }
     if (!is.null(dim(ali.mat)[1])) {
          if (dim(ali.mat)[1] != 0) {
               if (k <= dim(ali.mat)[1]) {
                    SP = ali.mat[order(-TM)[1:k], ]
                    SP = SP[!is.na(SP[, 1]), ]
               } else {
                    SP = ali.mat
               }
               if (dim(SP)[1] < k) {
                    print('Less Starting points selected')
               }
          } else {
               SP = matrix(0, 1, 1)
          }
     } else {
          SP = ali.mat
          print('Less Starting points selected')
     }
     return(SP)
}
