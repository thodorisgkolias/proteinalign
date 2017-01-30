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



gpa_C<-function(x,rescale=1,reflect=0,tol1=1e-5,tol2=tol1) {
     #input is the 3 dimensional matrix as in procGPA, tolerances are the same
     #rescale and reflect are obvious if values are 1 then the same as TRUE otherwise FALSE
     
     z=list()
     d<-dim(x)
     m<-d[2]
     k<-d[1]
     n<-d[3]
     x0<-x
     #x0<-aperm(x,c(2,1,3));
     a<-.C("ctorrgpa",as.double(x0),as.integer(k),as.integer(m),
           as.integer(n),as.integer(rescale),as.integer(reflect),
           as.double(tol1),as.double(tol2),out=x0)$out
     rm(x0);
     #z$rot<-aperm(a1,c(2,1,3));
     z$rot<-array(a,c(k,m,n))
     rm(a);
     z$mshape<-apply(z$rot,c(1,2),mean);
     
     return(z);
}


diam_C<-function(x,preshape=0,reflect=0) {
     #this produces the sample diameter, 
     #in particular gives the riemdist between two configurations if dim[3]=2
     # x is a 3-d matrix
     z=list()
     d<-dim(x)
     m<-d[2]
     k<-d[1]
     n<-d[3]
     #x0<-aperm(x,c(2,1,3));
     a<-.C("diam",as.double(x),as.integer(k),as.integer(m),as.integer(n),
           as.integer(preshape),as.integer(reflect),out=double(length(1)))$out
     #rm(x0);
     return(a);
}



riemdist_C<-function(x,y,preshape=0,reflect=0) {
     #riemdist between configurations x,y. If preshape=1 x,y are preshapes
     #reflect=1 then refleltion information is removed.

     Y<-array(0,c(dim(x)[1],dim(x)[2],2))
     Y[,,1]<-x
     Y[,,2]<-y
     s<-diam_C(Y,preshape,reflect)
     s
}

rgpa_C<-function(x,reflect=0,tol1=1e-10) {
     #input is the 3 dimensional matrix as in rgpa, tolerance is the same
     #rescale and reflect are obvious if values are 1 then the same as TRUE otherwise FALSE
     

     z=list()
     d<-dim(x)
     m<-d[2]
     k<-d[1]
     n<-d[3]
     #x0<-aperm(x,c(2,1,3));
     x0<-x
     a<-.C("rpackagergpa",as.double(x0),as.integer(k),as.integer(m),
           as.integer(n),as.integer(reflect),as.double(tol1),out=x0)$out
     rm(x0)
     #a1<-array(a,c(m,k,n));
     #z$rot<-aperm(a1,c(2,1,3));
     z$rot<-array(a,c(k,m,n))
     #z$mshape<-apply(z$rot,c(1,2),mean);
     rm(a)
     return(z)
}


sgpa_C<-function(x) {
     #input is the 3 dimensional matrix as in sgpa

     d<-dim(x)
     m<-d[2]
     k<-d[1]
     n<-d[3]
     #x0<-aperm(x,c(2,1,3));
     x0<-x
     a<-.C("rpackagesgpa",as.double(x),as.integer(k),as.integer(m),
           as.integer(n),out=x0)$out
     #rm(x0);
     #a1<-array(a,c(m,k,n));z<-aperm(a1,c(2,1,3));
     z<-array(a,c(k,m,n))
     #z$mshape<-apply(z$rot,c(1,2),mean);
     rm(a)
     return(z)
}

preshapetoicon_C<-function(x) {
     #input is the 2 dimensional matrix 

     d<-dim(x)
     m<-d[2]
     k<-d[1]
     x0<-aperm(x,c(2,1))
     a<-.C("preshapetoiconctor",as.double(x0),as.integer(k),as.integer(m),
           out=double((length(x)+m)))$out
     rm(x0);
     a1<-array(a,c(m,k+1))
     z<-aperm(a1,c(2,1))
     #z$mshape<-apply(z$rot,c(1,2),mean);
     rm(a)
     return(z)
}

cnt3_C<-function(x) {
     #input is the 3 dimensional matrix as in cnt3

     z=list()
     d<-dim(x)
     m<-d[2]
     k<-d[1]
     n<-d[3]
     #x0<-aperm(x,c(2,1,3));
     x0<-x
     a<-.C("rpackagercnt3",as.double(x),as.integer(k),as.integer(m),
           as.integer(n),out=x0)$out
     rm(x0)
     #z<-array(a,c(m,k,n));rm(a);z<-aperm(z,c(2,1,3));
     z<-array(a,c(k,m,n))
     #z$mshape<-apply(z$rot,c(1,2),mean);
     rm(a)
     return(z)
}


preshape_C<-function(x,rescale=1) {
     #input is the 2 dimensional matrix 

     d<-dim(x)
     m<-d[2]
     k<-d[1]
     x0<-aperm(x,c(2,1))
     a<-.C("preshapector",as.double(x0),as.integer(k),as.integer(m),
           as.integer(rescale),out=double((length(x)-m)))$out
     rm(x0)
     a1<-array(a,c(m,k-1))
     z<-aperm(a1,c(2,1))
     #z$mshape<-apply(z$rot,c(1,2),mean);
     rm(a)
     return(z)
}


###

vectorhelmert_C<-function(X){
     #function produces the matrix of preshapes
     #X is kxmxn output is Y a (k-1)xmxn matrix

     d<-dim(X)
     d[1]<-d[1]-1
     Y<-array(0,d)
     for (i in 1:d[3])
          Y[,,i]<-preshape_C(X[,,i])
     Y
}



