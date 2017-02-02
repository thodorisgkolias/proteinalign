
RMSD = function(X, Y) {
     ## computes the rmsd of two matrices
     n = dim(X)[1]
     res = sqrt(sum(Norm(X - Y)^2) / n )
     return(res)
}


SimilarityMetrics = function(data1, data2, sol) {
     L = dim(sol)[1]
     data_proc = array(NA, dim = c(dim(sol)[1], 3, 2))
     data_proc[, , 1] = data1[sol[, 1], ]
     data_proc[, , 2] = data2[sol[, 2], ]
     rot = GpaC(data_proc)$rot
     ## RMSD
     rmsd = round(RMSD(rot[, , 1], rot[, , 2]), 1)
     ## TMscore
     Ltar = dim(data1)[1]
     Lali = dim(data2)[1]
     di = apply(rot[, , 1] - rot[, , 2], 1, Norm)
     d01 = 1.24 * (Ltar - 15) ^ (1 / 3) - 1.8
     d02 = 1.24 * (Lali - 15) ^ (1 / 3) - 1.8
     TM1 = (1 / Ltar) * sum(1 / (1 + (di / d01) ^ 2))
     TM2 = (1 / Lali) * sum(1 / (1 + (di / d02) ^ 2))
     TMscore = c(TM1, TM2)
     ## SO
     m = min(Ltar, Lali)
     SO = round(100 * (sum(di < 3.5) / m), 2)
     #SPscore
     #dinew = di[which(di<=7)]
     #SPscore = (1/L^0.7)*sum((1/(1 + (dinew/3.5)^2)) - 0.2)
     
     return(list(rmsd = rmsd, TMscore = TMscore, SO = SO))
}


SimilarityTransf = function(X,Xrot) {
     n = dim(X)[3]
     m = dim(X)[2]
     R = array(, dim=c(m, m, n))
     location = matrix(, m, n)
     X_cnt = Cnt3C(X)
     Xrot_cnt = Cnt3C(Xrot)
     for (i in 1:n) {
          location[, i] = apply(X[, , i] - Xrot[, , i], 2, mean)
          SS = SvdNew(t(X_cnt[, , i]) %*% Xrot_cnt[, , i])
          R[, , i] = t(SS$v %*% t(SS$u))
     }
     out = list(location = location , R = R)
     return(out)
}
