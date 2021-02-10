a3gr <- function(A){
    # Givens decomosition A to Q R where R is in original A
    m <- dim(A)[1]
    n <- dim(A)[2]
    Q <- diag(m)  
    tol <- .Machine$double.eps^(1.5)
#  ?? should we use 2 or 1.5 or ??
    mn <- min(m, n)
#  STEP 1
    if (m == 1) {
      sol <- list(Q=Q, Rx=A)
      return(sol)
    }
    for (j in 1:mn){
       if (j+1 > m) break # end loop
#  C  STEP 2
      for (k in (j+1):m){
         # C  STEP 3
         C <- A[j,j]
         S <- A[k,j]
         B <- abs(C)
         if (abs(S) > B) {B <- abs(S)}
         if (B == 0.0) break # goto 90
         C <- C/B
         S <- S/B
         P <- sqrt(C*C+S*S)
#     C  STEP 4
         S <- S/P
#     C  STEP 5
         if (abs(S) <= tol) {break} # GOTO 90 Note: <= rather than <
#    C  STEP 6
         C <- C/P
#    C  STEP 7
         Pv <- A[j,]
         A[j,] <- C*Pv+S*A[k,]
         A[k,] <- -S*Pv+C*A[k,]
# C  STEP 8
         Pv <- Q[,j]
         Q[,j] <- C*Pv+S*Q[,k]
         Q[,k] <- -S*Pv+C*Q[,k]
# C  STEP 9   90    CONTINUE
      }
# C  STEP 10   100  CONTINUE
    }
    sol <- list(Q=Q, Rx=A)
    sol 
}
## C  TEST ALGORITHM 3
# dro3ca.R -- attempt to use implicit looping
# m <- as.numeric(readline("no. of rows="))
# n <- as.numeric(readline("no of columns="))
m<-5
n<-3
# build frank matrix
cat("Frank matrix A ",m," by ",n,"\n")
A <- matrix(0.0, nrow=m, ncol=n)
for (i in 1:m){
    for (j in 1:n){
        A[i,j] <- min(i,j)-2.0
    }
    if (i <= j) {A[i,i] <- i}
}
print(A)
Acopy <- A
grA <- a3gr(A)
Q <- grA$Q
cat("Q\n")
print(Q)
R <- grA$Rx
cat("R\n")
print(R)
test <- Q %*% R
cat("error =",max(abs(test-Acopy)),"\n")
