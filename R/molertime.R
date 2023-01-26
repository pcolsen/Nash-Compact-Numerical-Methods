molerfast <- function(n) {
# A fast version of `molermat'
A <- matrix(0, nrow = n, ncol = n)
j <- 1:n
for (i in 1:n) {
A[i, 1:i] <- pmin(i, 1:i) - 2
}
A <- A + t(A)
diag(A) <- 1:n
A
}

molermat<-function(n){
   A<-matrix(NA, nrow=n, ncol=n)
   for (i in 1:n){
      for (j in 1:n) {
          if (i == j) A[i,i]<-i
          else A[i,j]<-min(i,j) - 2
      }
   }
   A
}

axmolerfast <- function(x) {
# A fast and memory-saving version of A%*%x  
# For Moler matrix
n <- length(x)
j <- 1:n
ax <- rep(0, n)
for (i in 1:n) {
term <- x * (pmin(i, j) - 2)
ax[i] <- sum(term[-i]) 
}
ax <- ax + j*x
ax
}

# library(microbenchmark)
library(bench)
for (i in 1:20){
n<-100*i
# treg<-system.time(A<-molermat(n))
# tfast<-system.time(AA<-molerfast(n))
# if (! identical(A, AA)) stop("bad matrices")
tmoler<- bench::mark( A<-molermat(n),
        AA<-molerfast(n))
if (! identical(A, AA)) stop("bad matrices")
cat("n =",n,"\n")
print(tmoler)
cat("\n")
}

