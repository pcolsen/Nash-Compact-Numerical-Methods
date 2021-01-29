tstsvdR <- function(A, ssol){
   # J C Nash 20190326
   # test proposed SVD 
   # reconstruction
   U <- ssol$u
   V <- ssol$v
   S <- ssol$d
   AR <- (U %*% diag(S)) %*% t(V)
   ra <- max(abs(AR-A)) # absolute
   rr <- ra/max(abs(A)) # one relative measure
   nn <- dim(A)[2]
   mm <- dim(A)[1]
   # U orthog
   UTU <- t(U) %*% U
   UUT <- U %*% t(U) 
   # V orthog
   VTV <- t(V) %*% V
   VVT <- V %*% t(V) 
   autu <- max(abs(UTU-diag(nn)))
   auut <- max(abs(UUT-diag(mm)))
   avtv <- max(abs(VTV-diag(nn)))
   avvt <- max(abs(VVT-diag(nn)))
   res <- list(ra=ra, rr=rr, autu=autu, auut=auut, avtv=avtv, avvt=avvt)
   res
}