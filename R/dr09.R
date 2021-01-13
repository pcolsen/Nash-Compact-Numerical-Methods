A9 <- function(a, n){
    x <- rep(0, n)
    for (k in n:1){
        s=a[1]
        if (s <= 0){
          stop("A9: matrix is singular")
        }
       m<-1
       for (i in 2:n){
         q<-m; m<-m+i; t<-a[q+1]; x[i]<--t/s
         if (i > k) { x[i] <- -x[i]}
         for (j in (q+2):m){
           a[j-i]<-a[j]+t*x[j-q]
         }
       }
     q<-q-1; a[m]=1/s
     for (i in 2:n){a[q+i] <- x[i]}  
#    cat("iteration k:")
#     print(a)
    }
    a
}

FrankMat <- function(n){
  Amat <- matrix(0, nrow=n, ncol=n)
  for (i in 1:n){
     for (j in 1:i){
         Amat[i,j]<-j; Amat[j,i]<-j
     }
  }
    Amat
}

smat2vec <- function(Amat){
   n<-dim(Amat)[1]
   n2<-(n*(n+1)/2)
   svec = rep(0, n2)
   k <- 0
  for (i in 1:n){
    for (j in 1:i){
       k<-k+1
       svec[k]<-Amat[i,j]
    }
  }
  svec
}

svec2mat <- function(svec){
  n2<-length(svec)
  n <- (-1+sqrt(1+8*n2))/2
  Amat <- matrix(0, nrow=n, ncol=n)
  k <- 0
  for (i in 1:n){
    for (j in 1:i){
      k<-k+1
      Amat[j,i]<-Amat[i,j]<-svec[k]
    }
  }
  Amat
}

AA <- FrankMat(4)
vv <- smat2vec(AA)
vv
vinv<-A9(vv, 4)
vinv
print(vinv)
Ainv<-svec2mat(vinv)
print(Ainv)
print(Ainv %*% AA)
