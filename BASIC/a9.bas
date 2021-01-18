10 PRINT "ALGORITHM 9 - BAUER REINSCH INVERSION TEST"
20 N=100
40 DIM A(N*(N+1)/2),X(N)
45 LET N=4
50 GOSUB 1500
51 REM BUILD MATRIX IN A
60 GOSUB 1400 
61 REM PRINT IT
70 GOSUB 1000
71 REM INVERT
80 GOSUB 1400
81  REM PRINT
110 STOP
1000 REM ALG. 9 BAUER REINSCH INVERSION
1010 FOR K=N TO 1 STEP -1
1011     REM STEP 1
1020   S=A(1)
1021     REM STEP 2
1030   IF S<=0 THEN EXIT 1160
1031     REM STEP 3
1040   M=1
1041     REM STEP 4
1050   FOR I=2 TO N
1051     REM STEP 5
1060      Q=M
1061     M=M+I
1062     T=A(Q+1)
1063     X(I)=-T/S
1064     REM STEP 6
1070      IF I>K THEN X(I)=-X(I)
1071     REM STEP 7
1080      FOR J=Q+2 TO M
1081     REM STEP 8
1090         A(J-I)=A(J)+T*X(J-Q)
1100      NEXT J
1110   NEXT I
1111     REM STEP 9
1120   Q=Q-1
1121     A(M)=1/S
1122     REM STEP 10
1130   FOR I=2 TO N
1131     A(Q+I)=X(I)
1132     NEXT I
1133     REM STEP 11
1140 NEXT K
1141     REM STEP 12
1150 RETURN
1160 PRINT "MATRIX COMPUTATIONALLY INDEFINITE"
1170 STOP
1171     REM END ALG. 9
1400 PRINT "MATRIX A"
1410 FOR I=1 TO N
1420 FOR J=1 TO I
1430 PRINT A(I*(I-1)/2+J);
1440 NEXT J
1450 PRINT
1460 NEXT I
1470 RETURN
1500 REM FRANK MATRIX
1510 FOR I=1 TO N
1520 FOR J=1 TO I
1530 LET A(I*(I-1)/2+J)=J
1540 NEXT J
1550 NEXT I
1560 RETURN
