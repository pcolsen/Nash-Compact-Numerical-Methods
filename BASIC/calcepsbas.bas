100 PRINT "calceps.bas -- machine precision etc."
110 GOSUB 7000
120 PRINT  "B9 = ";B9;"  -- a large number"
130 PRINT "E5 = ";E5;"  -- a number for relative equality tests"
140 PRINT  "E9 = ";E9;"  -- the machine precision*1e+15,  E9=MIN (X : 1+X>1)"
150 PRINT "E6 = ";E6;"  -- the radix of arithmetic"
160 PRINT "J1 = ";J1;"  -- the number of radix digits in mantissa of FP numbers"
170 PRINT "E6^(-J1+1) =";(E6^(-J1+1))
180 SYSTEM
7000 REM ENVRON -- SUBROUTINE TO DETERMINE MACHINE PRECISION
7010 REM                AND SET SOME STANDARD VARIABLES
7020 REM INPUTS:  NONE
7030 REM
7040 REM OUTPUTS:
7050 REM    B9 -- a large number (currently 1E+35)
7060 REM    E5 -- a number for relative equality tests (currently 10)
7070 REM    E9 -- the machine precision,  E9=MIN (X : 1+X>1) 
7080 REM    E6 -- the radix of arithmetic
7090 REM    J1 -- the number of radis digits in mantissa of
7100 REM          floating-point numbers
7110 REM
7120 LET D1=1: REM use a variable to avoid constants when double
7130 REM       precision is invoked
7140 LET E5=10: REM arbitrary scaling for additive equality tests
7150 LET B9=1E+35: REM big number, not necessarily biggest possible
7160 LET E6=1: REM initial value for radix 
7170 LET E9=1: REM initial value for machine precision
7180 LET E9=E9/2: REM start of loop to decrease estimated machine precision
7190 LET D0=E6+E9: REM force storage of sum into a floating-point scalar
7200 IF D0>E6 THEN 7180: REM repeat reduction while (1+E9) > 1
7210 LET E9=E9*2: REM restore smallest E9 which gives (1+E9) > 1
7220 LET E6=E6+1: REM try different radix values
7230 LET D0=E6+E9
7240 IF D0>E6 THEN 7220: REM until a shift is observed
7250 LET J1=1: REM initial count of radix digits in mantissa
7260 LET E9=1: REM use radix for exact machine precision
7270 LET E9=E9/E6: REM loop while dividing by radix
7280 LET J1=J1+1: REM increment counter
7290 LET D0=D1+E9: REM add tp 1
7300 IF D0>D1 THEN 7270: REM test and repeat until equality
7310 LET E9=E9*E6: REM recover last value of machine precision
7320 LET J1=J1-1: REM and adjust the number of digits
7330 RETURN
