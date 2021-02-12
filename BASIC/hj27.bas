40 DIM T(25): REM DIMension only needed for POSTGEN, not HJ
1000 PRINT "Hooke and Jeeves -- 19851018, 19880809"
1004 PRINT #3,"Hooke and Jeeves -- 19851018, 19880809"
1008 REM CALLS:
1012 REM     FUNCTION F(B)  -- line 2000
1016 REM     ENVRON (computing  environment) -- line 7120
1020 REM
1024 REM INPUTS TO THE ROUTINE:
1028 REM    B() -- a vector of initial parameter estimates
1032 REM    S1  -- initial stepsize for search (set to 1 if zero)
1036 REM    S2  -- stepsize reduction factor, which is applied to
1040 REM           S1 when axial search fails (set to 0.1 if zero)
1044 REM    N   -- the number of parameters in the function F(B)
1048 REM    O( , ) -- masks and bounds 
1052 REM    I5  -- the number of masked parameters. must be zero if
1056 REM           masks are not used
1060 REM OUTPUT FROM THE ROUTINE:
1064 REM    X() -- a vector of final parameter estimates for the
1068 REM           values of the parameters which minimize the function.
1072 REM    F0  -- value of the function at the minimum
1076 REM    I8  -- number of gradient evaluations (unchanged)
1080 REM    I9  -- number of function evaluations
1084 REM 
1088 IF S1<=0 THEN LET S1=1: REM !! warning -- is variable undefined?
1092 IF S2<=0 THEN LET S2=.1: REM !!  ditto
1096 LET J8=1: REM no of fn eval before parameter display - mod 19880809
1100 LET J7=1: REM counter for parameter display
1104 GOSUB 7120: REM computing environment
1108 GOSUB 1440: REM copy B() into X() (lowest point so far)
1112 PRINT "STEP-SIZE =";S1
1116 PRINT "STEP-SIZE REDUCTION FACTOR =";S2
1120 PRINT #3,"STEP-SIZE =";S1
1124 PRINT #3,"STEP-SIZE REDUCTION FACTOR =";S2
1128 GOSUB 1456: REM compute function in F, set I3<>0 if not possible.
1132 IF I3<>0 THEN 1412
1136 PRINT "INITIAL FUNCTION VALUE =";F
1140 PRINT #3,"INITIAL FUNCTION VALUE =";F
1144 LET F1=F: REM store function value at base point
1148 LET F0=F: REM store lowest function value so far
1152 GOSUB 1252: REM axial exploratory search
1156 IF I6=2*(N-I5) THEN 1176: REM parameters unchanged in axial search
1160 IF F0>=F1 THEN 1176: REM test for a lower function value
1164 LET F1=F0: REM update function value at base
1168 GOSUB 1368: REM pattern move
1172 GOTO 1152: REM repeat axial search
1176 FOR J=1 TO N: REM is B() still the current base point?
1180 IF B(J)<>X(J) THEN 1192: REM test for changes in parameters
1184 NEXT J: REM in above test look for equality since B:=X at base
1188 GOTO 1208: REM reduce step-size as search has not reduced function
1192 GOSUB 1424: REM copy X into B (B() is now at the base point)
1196 PRINT " RETURN TO BASE POINT ";
1200 PRINT #3," RETURN TO BASE POINT ";
1204 GOTO 1152: REM try another axial search
1208 LET S1=S1*S2: REM reduce step-size
1212 PRINT
1216 PRINT I9,F0,"STEPSIZE=";S1
1220 PRINT #3,
1224 PRINT #3,I9,F0,"STEPSIZE=";S1
1228 GOSUB 1476: REM display parameters
1232 IF I6<2*(N-I5) THEN 1152: REM convergence test (no altered param's)
1236 PRINT
1240 PRINT #3,
1244 RETURN: REM function minimization complete
1248 REM axial exploratory search subroutine
1252 LET I6=0: REM counter for number of unchanged parameters
1256 FOR J=1 TO N
1260 IF O(J,3)=0 THEN 1344: REM parameter J is masked
1264 LET S3=B(J): REM store parameter value
1268 IF S3=O(J,2) THEN 1284: REM avoid change if at upper bound
1272 LET B(J)=S3+S1: REM step forward
1276 IF B(J)>O(J,2) THEN LET B(J)=O(J,2): REM bound -- assume +VE step
1280 IF B(J)+E5<>S3+E5 THEN 1292: REM test equality relative to E5
1284 LET I6=I6+1
1288 GOTO 1300: REM now try negative step
1292 GOSUB 1456: REM function evaluation
1296 IF F<F0 THEN 1340: REM test if function value < current lowest value
1300 IF S3=O(J,1) THEN 1328: REM no change if at lower bound
1304 LET B(J)=S3-S1: REM step backward
1308 IF B(J)<O(J,1) THEN LET B(J)=O(J,1): REM lower bound
1312 IF B(J)+E5=S3+E5 THEN 1328: REM test equality
1316 GOSUB 1456: REM function evaluation
1320 IF F<F0 THEN 1340
1324 GOTO 1332
1328 LET I6=I6+1: REM count number of times parameter unchanged
1332 LET B(J)=S3: REM restore original parameter (not by addition!!)
1336 GOTO 1344
1340 LET F0=F: REM store new lowest function value
1344 NEXT J
1348 PRINT " AXIAL SEARCH F0=";F0 
1352 PRINT #3," AXIAL SEARCH F0=";F0 
1356 GOSUB 1476: REM print parameters
1360 RETURN: REM end axial search
1364 REM PATTERN MOVE
1368 FOR J=1 TO N
1372 IF O(J,3)=0 THEN 1396: REM masked parameter
1376 LET S3=2*B(J)-X(J): REM element of new base point
1380 LET X(J)=B(J): REM store current point
1384 IF S3<O(J,1) THEN LET S3=O(J,1): REM lower bound check
1388 IF S3>O(J,2) THEN LET S3=O(J,2): REM upper bound check
1392 LET B(J)=S3: REM store new base point
1396 NEXT J
1400 PRINT " PMOVE ";
1404 PRINT #3," PMOVE ";
1408 RETURN: REM end pattern move
1412 PRINT "FUNCTION NOT COMPUTABLE AT INITIAL POINT"
1416 PRINT #3,"FUNCTION NOT COMPUTABLE AT INITIAL POINT"
1420 STOP
1424 FOR J=1 TO N: REM copy X into B
1428 LET B(J)=X(J)
1432 NEXT J
1436 RETURN
1440 FOR J=1 TO N: REM copy B into X
1444 LET X(J)=B(J)
1448 NEXT J
1452 RETURN
1456 LET I3=0: REM compute function -- reset failure flag
1460 GOSUB 2000: REM user routine
1464 IF I3<>0 THEN LET F=B9: REM large value assigned
1468 LET I9=I9+1: REM function evaluation counter
1472 RETURN
1476 IF J8=0 THEN RETURN: REM no parameter display
1480 IF I9<J7*J8 THEN RETURN: REM check if to be printed
1484 LET J7=INT(I9/J8)+1: REM parameter display control
1488 PRINT "parameters";
1492 PRINT #3,"parameters";
1496 FOR J=1 TO N
1500 LET Q$=""
1504 IF B(J)-O(J,1)<=E9*(ABS(O(J,1))+E9) THEN LET Q$="L": REM lower bd
1508 IF O(J,2)-B(J)<=E9*(ABS(O(J,2))+E9) THEN LET Q$="U": REM upper bd
1512 IF O(J,3)=0 THEN LET Q$="M": REM masked 
1516 PRINT "  ";B(J);Q$;
1520 PRINT #3,"  ";B(J);Q$;
1524 IF 5*INT(J/5)<>J THEN 1544
1528 PRINT
1532 PRINT "          ";
1536 PRINT #3,
1540 PRINT #3,"          ";
1544 NEXT J
1548 PRINT
1552 PRINT #3,
1556 RETURN
