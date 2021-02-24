cat("R timer for exp(sin(cos)))\n")
# Rtimetest.R
n <- 2345678
timeexsico <- function(n){
  y <- 0.0
  for (i in 1:n){
    x<-exp(sin(cos(i)))
    y<-y+x
  }
  y
}

system.time(timeexsico(n))
library(microbenchmark)
microbenchmark(timeexsico(n))
cat("Done!\n")