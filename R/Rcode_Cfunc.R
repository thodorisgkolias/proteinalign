## rotation integral
myrotsim.laplaceC = function(l) {
     res = .Call("dointeg2", as.double(l), out = as.double(rep(0, 4)))$out
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
     a = .Call("EM_C",as.double(x0), as.integer(k), as.integer(m), as.integer(n),
            out = x0, sig = 0,lik = 0)
     rm(x0)
     a$out = a$out[-k, , ]
     z$rot = array(a$out, c(k-1, m, n))
     z$mshape = apply(z$rot, c(1, 2), mean)
     return(list(mean = t(z$mshape), sigma = a$sig, lik = a$lik,
                 rot.data = aperm(z$rot, c(2, 1, 3)) ))
}
