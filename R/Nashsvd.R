Nashsvd <- function(A, MaxRank=0, cyclelimit=25, trace = 0, rotnchk=0.3) {
## Nashsvd.R -- An attempt to remove tolerances from Nash & Shlien algorithm 190327
# Partial svd by the one-sided Jacobi method of Nash & Shlien
#  Computer Journal 1987  30(3), 268-275
#  Computer Journal 1975  18(1) 74-76
  if (cyclelimit < 6) {
      warning("Nashsvd: You cannot set cyclelimit < 6 without modifying the code")
      cyclelimit <- 6 # safety in case user tries smaller
  }
  m <- dim(A)[1]
  n <- dim(A)[2]
  if (MaxRank <= 0) MaxRank <- n
  EstColRank <- n # estimated column rank 
  # Note that we may simply run algorithm to completion, or fix the
  # number of columns by EstColRank. Need ?? to fix EstColRank=0 case.??
  V <- diag(nrow=n) # identity matrix in V
  if (is.null(EstColRank)) {EstColRank <- n } # Safety check on number of svs
  z <- rep(NA, n) # column norm squares -- safety setting
  keepgoing <- TRUE
  SweepCount <- 0
  while (keepgoing) { # main loop of repeating cycles of Jacobi
    RotCount <- 0
    SweepCount <- SweepCount + 1 
    if (trace > 1) cat("Sweep:", SweepCount,"\n")
##    if (EstColRank == n) { EstColRank <- n - 1 } # safety
    for (jj in 1:(EstColRank-1)) { # left column indicator
       for (kk in (jj+1): n) { # right hand column
         p <- q <- r <- 0.0 # 
         oldjj <- A[,jj]
         oldkk <- A[,kk]
         p <- as.numeric(crossprod(A[,jj], A[,kk]))
         q <- as.numeric(crossprod(A[,jj], A[,jj]))
         r <- as.numeric(crossprod(A[,kk], A[,kk]))
         if (trace > 2) cat(jj," ",kk,": pqr",p," ",q," ",r," ")
         z[jj]<-q
         z[kk]<-r
         if (q >= r) { # in order, so can do test of "convergence" -- change to 0.2 * abs(p) for odd cases
            if ( (as.double(z[1]+q) > as.double(z[1]) ) && (as.double(rotnchk*abs(p)+q) > as.double(q)) )  {
              RotCount <- RotCount + 1
              p <- p/q
              r <- 1 - (r/q)
              vt <- sqrt(4*p*p +r*r)
              c0 <- sqrt(0.5*(1+r/vt))
              s0 <- p/(vt*c0)
              # rotate
              cj <- A[,jj]
              ck <- A[,kk]
              A[,jj] <- c0*cj + s0*ck
              A[,kk] <- -s0*cj + c0*ck
              cj <- V[,jj]
              ck <- V[,kk]
              V[,jj] <- c0*cj + s0*ck
              V[,kk] <- -s0*cj + c0*ck
            } else {
              if (trace > 2) cat(" NO rotn ") 
            }
         } else { # out of order, must rotate
            if (trace > 2) cat("|order|")
            RotCount <- RotCount + 1
            p <- p/r
            q <- (q/r) - 1.0
            vt <- sqrt(4*p*p +q*q)
            s0 <- sqrt(0.5*(1-q/vt))
            if (p < 0) { s0 <- -s0 }
            c0 <- p/(vt*s0)
            # rotate
            cj <- A[,jj]
            ck <- A[,kk]
            A[,jj] <- c0*cj + s0*ck
            A[,kk] <- -s0*cj + c0*ck
            cj <- V[,jj]
            ck <- V[,kk]
            V[,jj] <- c0*cj + s0*ck
            V[,kk] <- -s0*cj + c0*ck
         } # end q >= r test
         nup <- as.numeric(crossprod(A[,jj], A[,kk]))
#         nuq <- as.numeric(crossprod(A[,jj], A[,jj]))
#         nur <- as.numeric(crossprod(A[,kk], A[,kk]))
         if (trace > 2) cat("    new: p= ",nup," Rel:",nup*nup/z[1],"\n")
       } # end kk
    } # end jj
    if (trace > 0) {cat("End sweep ", SweepCount,"  No. rotations =",RotCount,"\n")}
    if (trace > 2) tmp <- readline("cont.?")
    while ( (EstColRank >= 3) && (as.double(sqrt(z[EstColRank])+sqrt(z[1]) == as.double(sqrt(z[1])) ))) {
    # ?? Why can we not use 2? Or do we need at least 2 cols
        EstColRank <- EstColRank - 1
        if (trace > 0) {cat("Reducing rank to ", EstColRank,"\n")} # ?? can do this more cleanly
    } # end while for rank estimation
    ## Here may want to adjust for MaxRank. How??
    if (MaxRank < EstColRank) {
       if (trace > 0) {
        cat("current estimate of sv[",MaxRank,"/sv[1] =",sqrt(z[MaxRank]/z[1]),"\n")
        cat("reducing rank by 1\n")
       }
       EstColRank <- EstColRank - 1
    }
    if ( SweepCount >= cyclelimit) { 
         if (trace > 0) cat("Cycle limit reached\n")
         keepgoing <- FALSE
    } 
    if (RotCount == 0) {
        if (trace > 1) cat("Zero rotations in cycle\n")
        keepgoing <- FALSE
    }
  } # End main cycle loop
  z <- sqrt(z)
  A <- A %*% diag(1/z)
  ans <- list( d = z, u = A, v=V, cycles=SweepCount, rotations=RotCount)
  ans
} # end partsvd()
