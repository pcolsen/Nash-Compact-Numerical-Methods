cat("calceps.R -- machine precision etc.\n")
envron<-function(){
  
# 7000 REM  ENVRON -- SUBROUTINE TO DETERMINE MACHINE PRECISION
# 7010 REM                AND SET SOME STANDARD VARIABLES
# 7020 REM INPUTS:  NONE
# 7030 REM
# 7040 REM OUTPUTS:
#   7050 REM    B9 -- a large number (currently 1E+35)
# 7060 REM    E5 -- a number for relative equality tests (currently 10)
# 7070 REM    E9 -- the machine precision,  E9=MIN (X : 1+X>1) 
# 7080 REM    E6 -- the radix of arithmetic
# 7090 REM    J1 -- the number of radis digits in mantissa of
# 7100 REM          floating-point numbers
# 7110 REM
D1<-as.numeric(1) # use a variable to avoid constants when double
#      precision is invoked
E5 <- 10 #arbitrary scaling for additive equality tests
B9 <- 1E+35 #big number, not necessarily biggest possible
cat("B9 = ",B9,"  -- a large number","\n")
cat("E5 = ",E5,"  -- a number for relative equality tests","\n")
E6 <- 1 #initial value for radix 
E9 <- 1 #initial value for machine precision
D0 <- E6
repeat{
  E9 <- E9/2 #start of loop to decrease estimated machine precision
  D0 <- E6+E9 #force storage of sum into a floating-point scalar
#  cat("E9temp=",E9,"\n")
#  readline("cont.")
  if (D0 <= E6) break
} #repeat reduction while (1+E9) > 1
E9 <- E9*2 #restore smallest E9 which gives (1+E9) > 1
cat( "TRIAL EPS=",E9,"\n")
repeat {
   E6 <- E6+1 #try different radix values
   D0 <- E6+E9
   if (D0 <= E6) break
} #until a shift is observed
cat("ESTIMATED RADIX=",E6,"\n")
J1 <- 1 #initial count of radix digits in mantissa
E9 <- 1 #use radix for exact machine precision
repeat {  
  E9 <- E9/E6 #loop while dividing by radix
  J1 <- J1+1 #increment counter
  D0 <- D1+E9 #add tp 1
  if (D0<=D1) break #test and repeat until equality
}

E9 <- E9*E6 #recover last value of machine precision
J1 <- J1-1 #and adjust the number of digits
mpvals<-list(eps=E9, radix=E6, ndigits=J1)
}
mp<-envron()
E9<-mp$eps
E6<-mp$radix
J1<-mp$ndigits
cat("E9 = ",E9,"  -- the machine precision,  E9=MIN (X  1+X>1)\n")
cat( "E9*1E+16=",E9*1E+16,"\n")
cat( "E6 = ",E6,"  -- the radix of arithmetic","\n")
cat( "J1 = ",J1,"  -- the number of radix digits in mantissa of FP numbers","\n")
cat( "E6^(-J1+1) * 1E+16=",1E+16*(E6^(-J1+1)),"\n")
cat(".Machine$double.eps=",.Machine$double.eps,"\n")

