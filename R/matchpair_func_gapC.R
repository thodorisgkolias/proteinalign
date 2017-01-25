##############################################
#### Fucntions for the matching algorithm #### 
### Updated : 22/12/2016                   ###
##############################################

#### Read data ####

data_extract = function(data, atom = 'CA') {
     ## data must be an output from readpdb.
     ## selects the amino acid sequence and the coordinates columns.
     temp = data[data[, 3] == atom, c(4, 7, 8, 9)]
     return(temp)
}


split_str = function(st) {
     temp = unlist(strsplit(st, ' '))
     temp2 = temp[temp != '']
     return(temp2)	
}

amino_transf = function(data) {
     ## Transforms the names of the amino acid sequence.
     
     freqs1 = c('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
                'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
                'TYR', 'VAL')
     freqs2 = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                'F', 'P', 'S', 'T', 'W', 'Y', 'V')
     freqs = cbind(freqs1, freqs2)
     for (i in 1:length(levels(data$amino))) {
          aa = levels(data$amino)[i]
          rows = which(levels(data$amino)[i] == freqs[, 1])
          levels(data$amino) = sub(paste('^', aa, '$', sep =''),
                                   as.character(freqs[rows, 2]), 
                                   levels(data$amino))
     }
     return(data)
}

readpdb = function(pdb, atom = 'CA', chain = 'A') {
     if (nchar(pdb) != 4) {
          stop('No PDB file found')
     }
     ind1 = c(1, 7, 13, 17, 18, 22, 23, 27, 31, 39, 47, 55, 61, 73, 77, 79)
     ind2 = c(6, 11, 16, 17, 21, 22, 26, 27, 38, 46, 54, 60, 66, 76, 78, 80)
     file = paste('http://www.rcsb.org/pdb/files/', pdb, '.pdb', sep = '')
     data1 = readLines(file, n =-1)
     choose = substring(data1, 1, 6)
     atoms = data1[choose == 'ATOM  ']
     clean.data = lapply(atoms, substring, ind1, ind2)
     newd = lapply(clean.data, split_str)
     rm_col = which(lapply(newd, length) != 12)
     if(length(rm_col) != 0) {
          for (i in rm_col) {
               newd[[i]] = newd[[i]][-7]
          }
     }
     newmat = matrix(unlist(newd), ncol = 12, byrow = TRUE)
     newmat = newmat[newmat[, 3] == atom, ]
     if (newmat[1, 5] != 'A') {
          chain = newmat[1, 5]
     }
     newmat = newmat[newmat[, 5] == chain, ]
     colnames(newmat) = c('type', 'seq', 'residue', 'amino', 'chain', 'res_no',
                          'x', 'y', 'z', 'occup', 'temp', 'elem_name')
     out = data.frame(newmat)
     out[, 7] = as.numeric(levels(out[, 7]))[out[, 7]]
     out[, 8] = as.numeric(levels(out[, 8]))[out[, 8]]
     out[, 9] = as.numeric(levels(out[, 9]))[out[, 9]]
     cat('', 'Loading...', '\n')
     cat('', split_str(data1[choose == 'HEADER'])[-1], '\n')
     
     out = data_extract(data = amino_transf(data = out),  atom = atom)
     return(out)
}



readdata = function(dir){ 
     ## Read pdb format data from a directory
     
     ind1 = c(1, 7, 13, 17, 18, 22, 23, 27, 31, 39, 47, 55, 61, 73, 77, 79)
     ind2 = c(6, 11, 16, 17, 21, 22, 26, 27, 38, 46, 54, 60, 66, 76, 78, 80)
     data1 = readLines(dir, n =-1)
     choose = substring(data1, 1, 6)
     atoms = data1[choose == 'ATOM  ']
     clean.data = lapply(atoms, substring, ind1, ind2)
     newd = lapply(clean.data, split_str)
     ll = lapply(newd, length)
     ll = unlist(ll)
     for (jj in 1:length(ll)){
          newd[[jj]] = newd[[jj]][1:min(ll)]
     }
     newmat = matrix(unlist(newd), ncol = ll, byrow = TRUE)
     out = data.frame(newmat)
     out[, 6] = as.numeric(levels(out[, 6]))[out[, 6]]
     out[, 7] = as.numeric(levels(out[, 7]))[out[, 7]]
     out[, 8] = as.numeric(levels(out[, 8]))[out[, 8]]
     out = amino_transf(out[out[, 3] == 'CA',c(4, 6, 7, 8)])
     return(out)
}


pam.transf = function(pam) {
     res = round(10^(pam/10), 6)
     return(res)
}



#### Likelihood functions ####

riemd = function(X, Y) {
     ## Riemannian distance
     res = ssriemdist(X, Y, FALSE)
     return(res)
}

log.lik.nodD = function(X, mu, sigma, p, M, A) {
     ## Computes the size and shape likelihood of a matched matrix X
     ## from a mean mu inside a volume of a cube A.
     m = dim(X)[1]
     Q = (p-1) * m
     RDs = X[, 1:(p-1)] %*% t(mu[, 1:(p-1)]) / sigma^2
     ss = svd.m_new(RDs)
     ss$d = c(1, 1, sign(det(ss$v %*% t(ss$u)))) * ss$d
     lik = myrotsim.laplaceC(ss$d)[1]
     
     res = (p - M) * log(A) - (Q / 2) * log(2 * pi * sigma^2) - 
          (norm(X[, 1:(p-1)])^2 + norm(mu[, 1:(p-1)])^2) / (2 * sigma^2) + 
          lik
     return(res)
}



#### Selection of starting points ####

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
     res1 = msa(x0, method = 'Muscle',
                type = 'protein')
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
          proc = procGPA(data.proc)
          euc.dist = NULL
          for (i in 1:dim(ali.mat)[1]) {
               euc.dist[i] = dist(rbind(proc$rotated[i, , 1],
                                        proc$rotated[i, , 2]))
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

#### Rotation integration in C ####
## rotation integral
myrotsim.laplaceC = function(l) {
     res = .C("dointeg2", as.double(l), out = as.double(rep(0, 4)))$out
     return(res)
}
## EM 
EM_C<-function(x){
     # x is an array(k,m,n)
     z = list()
     d = dim(x)
     m = d[2]
     k = d[1]
     n = d[3]
     x0 = x
     a = .C("EM_C",as.double(x0), as.integer(k), as.integer(m), as.integer(n),
            out = x0, sig = 0,lik = 0)
     rm(x0)
     a$out = a$out[-k, , ]
     z$rot = array(a$out, c(k-1, m, n))
     z$mshape = apply(z$rot, c(1, 2), mean)
     return(list(mean = t(z$mshape), sigma = a$sig, lik = a$lik,
                 rot.data = aperm(z$rot, c(2, 1, 3)) ))
}

### Matching functions ####
norm = function(X) {
     Y = sqrt(sum(diag( X %*% t(X) ) ) )
     Y
}

svd.m_new = function(X) {
     SS = svd(X)
     SS$d = SS$d * c(1, 1, sign( det(SS$v %*% t(SS$u)) ) )
     SS$u = SS$u %*% diag(c(1, 1, sign(det(t(SS$u))))) 
     SS$v = SS$v %*% diag(c(1, 1, sign(det(t(SS$v)))))
     SS
}


indrows = function(comb2, mat) {
     ss1 = which(rownames(mat) == comb2[1])
     return(ss1)
}

indcols = function(comb2, mat) {
     ss1 = which(colnames(mat) == comb2[2])
     return(ss1)
}

vol.func = function(data) {
     vol = prod(apply(data, 2, max))
     return(vol)
}


PAM_log_lik = function (seq1, seq2, p, PAM) {
     ## As in fallaize
     score = 0
     if (!is.null(PAM)) {
          for (i in 1:p) {
               amino_row = which(seq1[i] == rownames(PAM))
               amino_col = which(seq2[i] == rownames(PAM))
               score = score + log(PAM[amino_row, amino_col])
          }
     }
     return(score)
}

lik.match.gapC = function(pair = 'NULL', data, seq, matched, vol, PAM, 
                          gap_open, gap_ext)	{
     ## Computes the likelihood of two matched molecules
     ## If pair != NULL the computes the likelihood of two 
     ## matched molecules with the extra pair included
     M1  = sapply(data,dim)[1, ]
     N = length(M1)
     m = dim(data[[1]])[2]
     hist = matched
     hist = rbind(hist, pair)
     pop = sapply(M1, seq_len)
     Xnew = list()
     p = dim(hist)[1]
     seq_new = list()
     for (i in 1:N) {
          matched.Xrows = hist[, i]
          unmatched.Xrows = which(!pop[[i]] %in% matched.Xrows)
          Xnew[[i]] = data[[i]][c(matched.Xrows, unmatched.Xrows), ]
          seq_new[[i]] = seq[[i]][c(matched.Xrows, unmatched.Xrows)]
          
     }
     EM.data = lapply(Xnew, t)
     EM.data2 =array(NA,dim=c(3, p, N))
     for (i in 1:N){
          EM.data2[, , i] = EM.data[[i]][, 1:p]
     }
     EM.data2 = aperm(EM.data2, c(2, 1, 3))
     res = EM_C(EM.data2)
     sigma = res$sigma
     lik = NULL
     for (i in 1:N){
          lik[i] = log.lik.nodD(res$rot.data[, , i], res$mean, sigma, p, 
                                M1[i], A = vol)
     }
     lik_str = sum(lik, na.rm = TRUE)
     lik_seq = PAM_log_lik(seq_new[[1]], seq_new[[2]], p, PAM)
     gap = gap_penalty(hist, gap_open, gap_ext, M1)
     lik0 = lik_str + lik_seq + gap
     return(lik0)	
}


hungarian.step = function(data, seq, lik, lik0, matched, rows2, vol, PAM, 
                          gap_open, gap_ext, comb2 = NULL) {
     ## Applies the Hungarian method for the matrix lik
     ## and returns the final alignment after testing if
     ## a pair can be added
     row.len = sapply(rows2, length)
     mat = matrix(lik, nrow = row.len[1], ncol = prod(row.len[-1]) )
     rownames(mat) = rows2[[1]]
     temp = expand.grid(rows2[-1])
     colnames(mat) = 1:prod(row.len[-1])
     if (all(mat < 0)){
          y = solve_LSAP(-mat)
     } else {
          y = solve_LSAP(mat, maximum = TRUE)
     }
     #     print(mem_used())
     ind = cbind(1:dim(mat)[1], y)
     ali = cbind(as.numeric(rownames(mat)), as.numeric(colnames(mat)[y]))
     liks = mat[ind]
     tt = cbind(ali, liks)
     tt = tt[order(-tt[, 3]), ]
     if (!is.numeric(dim(tt))) {
          pair = as.numeric(c(tt[1], temp[tt[2], ]))
          liknew = lik.match.gapC(pair, data, seq = seq, matched, vol,
                                  PAM = PAM, gap_open = gap_open, 
                                  gap_ext = gap_ext)
          if (liknew > lik0) {
               lik0 = liknew
               matched = rbind(matched, pair)
          }
     } else {
          k = 1
          liknew = NULL
          while (k <= dim(tt)[1]) {
               pair = c(tt[k, 1], as.numeric(temp[tt[k, 2], ]))
               liknew = lik.match.gapC(pair, data, seq = seq, matched, vol, 
                                       PAM = PAM, gap_open = gap_open,
                                       gap_ext = gap_ext)
               if (liknew > lik0) {
                    lik0 = liknew
                    matched = rbind(matched, pair)
                    k = k + 1
               } else {
                    break
               }
          }
     }
     return(list(lik0 = lik0, hist = matched))
}

remove.step = function(data, seq, matched, vol, lik0, PAM, gap_open, gap_ext) {
     remove = TRUE
     hist = matched
     while (remove) {
          lik.remove = NULL
          for (i in 1:dim(hist)[1]) {
               remove.pair = as.vector(hist[i, ])
               hist2 = hist[-i, ]
               lik.remove[i] = lik.match.gapC(pair = NULL, data = data, seq=seq, matched = hist2, 
                                           vol = vol, PAM = PAM, 
                                           gap_open = gap_open, 
                                           gap_ext = gap_ext)
          }
          if (max(lik.remove) > lik0) {
               print("Remove")
               pos.remove = which.max(lik.remove)
               lik0 = max(lik.remove)
               hist = hist[-pos.remove, ]
               pnew = dim(hist)[1]
               remove = TRUE
          } else {
               print("No more to remove")
               remove = FALSE
               pnew = dim(hist)[1]
          }
     }     
     
     return(list(lik0 = lik0, hist = hist))
}


gap_penalty = function(seqs, gap_open, gap_ext, len) {
     seq = seqs[order(seqs[, 1]), ]
     
     seq1 = c(0, seq[, 1], len[1] + 1)
     seq2 = c(0, seq[, 2], len[2] + 1)
     
     seq1 = diff(seq1)
     seq2 = diff(seq2)
     ff = abs(c(seq1, seq2))
     gap = 0
     for (i in 1:length(ff)) {
          if(ff[i] != 1) {
               gap = gap + gap_open + (ff[i] - 2) * gap_ext
          }
     }
     return(-gap)
}

hungarian.step1 = function(data, seq, lik, lik0, matched, comb2, rows2, vol,
                           PAM, gap_open, gap_ext) {
     ## Applies the Hungarian method for the matrix lik
     ## and returns the final alignment after testing if
     ## a pair can be added, when restrict search is TRUE
     
     row.len = sapply(rows2, length)
     if (all(lik < 0)){
          mat = matrix(-10^5, nrow = row.len[1], ncol = prod(row.len[-1]) )
     } else {
          mat = matrix(10^5, nrow = row.len[1], ncol = prod(row.len[-1]) )
          
     }   
     rownames(mat) = rows2[[1]]
     temp = expand.grid(rows2[-1])
     colnames(mat) = rows2[[2]]
     p1 = apply(comb2, 1, indrows, mat = mat)
     p2 = apply(comb2, 1, indcols, mat = mat)
     for (i in 1:dim(comb2)[1]) {
          mat[p1[i], p2[i]] = lik[i]
     }
     if (all(mat < 0)){
          y = solve_LSAP(-mat)
     } else {
          y = solve_LSAP(mat, maximum = TRUE)
     }
     ind = cbind(1:dim(mat)[1], y)
     ali = cbind(as.numeric(rownames(mat)), as.numeric(colnames(mat)[y]))
     liks = mat[ind]
     tt = cbind(ali, liks)
     tt = tt[order(-tt[, 3]), ]
     if (!is.numeric(dim(tt))) {
          pair = as.numeric(c(tt[1], temp[tt[2], ]))
          liknew = lik.match.gapC(pair, data, seq = seq, matched, vol, 
                                  PAM = PAM, gap_open = gap_open, 
                                  gap_ext = gap_ext)
          if (liknew > lik0) {
               lik0 = liknew
               matched = rbind(matched, pair)
          }
     } else {
          k = 1
          liknew = NULL
          while (k <= dim(tt)[1]) {
               pair = c(tt[k, 1], tt[k, 2])
               temp = rbind(matched, pair)
               temp = temp[order(temp[, 1]), ]
               if (all(diff(temp[, 1]) > 0)  & all(diff(temp[, 2]) > 0)) {
                    liknew = lik.match.gapC(pair, data, seq = seq, matched,
                                            vol, PAM = PAM, gap_open = gap_open,
                                            gap_ext = gap_ext)
                    if (liknew > lik0) {
                         lik0 = liknew
                         matched = rbind(matched, pair)
                         k = k + 1
                    } else {
                         break
                    }
               } else {
                    k = k + 1
               }
          }
     }
     return(list(lik0 = lik0, hist = matched))
}


#### Similarity metrics ####

RMSD = function(X, Y) {
     ## computes the rmsd of two matrices
     n = dim(X)[1]
     res = sqrt(sum(norm(X - Y)^2) / n )
     return(res)
}


Similarity.metrics = function(data1, data2, sol) {
     L = dim(sol)[1]
     data.proc = array(NA, dim = c(dim(sol)[1], 3, 2))
     data.proc[, , 1] = data1[sol[, 1], ]
     data.proc[, , 2] = data2[sol[, 2], ]
     rot = procGPA(data.proc)$rot
     ## RMSD
     rmsd = round(RMSD(rot[, , 1], rot[, , 2]), 1)
     ## TMscore
     Ltar = dim(data1)[1]
     Lali = dim(data2)[1]
     di = apply(rot[, , 1] - rot[, , 2], 1, norm)
     d01 = 1.24 * (Ltar - 15)^(1 / 3) - 1.8
     d02 = 1.24 * (Lali - 15)^(1 / 3) - 1.8
     TM1 = (1 / Ltar) * sum(1 / (1 + (di / d01)^2))
     TM2 = (1 / Lali) * sum(1 / (1 + (di / d02)^2))
     TMscore = c(TM1, TM2)
     ## SO
     m = min(Ltar, Lali)
     SO = round(100 * (sum(di < 3.5) / m), 2)
     #SPscore
     #dinew = di[which(di<=7)]
     #SPscore = (1/L^0.7)*sum((1/(1 + (dinew/3.5)^2)) - 0.2)
     
     return(list(rmsd = rmsd, TMscore = TMscore, SO = SO))
}

