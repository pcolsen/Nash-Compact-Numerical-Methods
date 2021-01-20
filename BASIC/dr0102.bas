5 PRINT "dr0102.bas -- Nashlib Alg 01 and 02 driver"
10 PRINT "from ENHSVA APR 7 80 -- MOD 850519, remod 210113"
20 LET E1=1.0E-7 
30 PRINT "ONE SIDED TRANSFORMATION METHOD FOR REGRESSIONS VIA"
40 PRINT "THE SINGULAR VALUE DECOMPOSITION -- J.C.NASH 1973,79"
150 LET M=4
160 LET N=3
210 DIM Y(M,N+1),A(M,N),T(N,N),G(N),X(N),Z(N),U(N),B(M)
220 DIM F$(10)
230 LET F$="K"
236 PRINT "Prep matrix and RHS"
240 LET Y(1,1)=5
241 LET Y(1,2)=1.0E-6
242 LET Y(1,3)=1
243 LET B(1)=1
250 LET Y(2,1)=6
251 LET Y(2,2)=0.999999
252 LET Y(2,3)=1
253 LET B(2)=2
260 LET Y(3,1)=7
261 LET Y(3,2)=2.00001
262 LET Y(3,3)=1
263 LET B(3)=3
270 LET Y(4,1)=8
271 LET Y(4,2)=2.9999
272 LET Y(4,3)=1
273 LET B(4)=4
500 FOR I=1 TO M
510 FOR J=1 TO N-1
520 LET A(I,J)=Y(I,J)
530 NEXT J
535 quit
540 LET A(I,N)=E3
550 NEXT I
560 LET E2=N*N*E1*E1
570 PRINT 
580 FOR I=1 TO N
590 FOR J=1 TO N
600 LET T(I,J)=0
610 NEXT J
620 LET T(I,I)=1
630 NEXT I
640 LET I9=0
650 IF N=1 THEN GOTO 1150
660 LET N2=N*(N-1)/2
670 LET N1=N-1
680 LET N9=N2
690 LET I9=I9+1
700 FOR J=1 TO N1
710 LET J1=J+1
720 FOR K=J1 TO N
730 LET P=0
740 LET Q=0
750 LET R=0
760 FOR I=1 TO M
770 LET P=P+A(I,J)*A(I,K)
780 LET Q=Q+A(I,J)*A(I,J)
790 LET R=R+A(I,K)*A(I,K)
800 NEXT I
810 IF Q>=R THEN GOTO 850
820 LET C=0
830 LET S=1
840 GOTO 920
850 IF (Q*R)<=0 THEN GOTO 1040
860 IF P*P/(Q*R)<E2 THEN GOTO 1040
870 LET Q=Q-R
880 LET P=2*P
890 LET V1=SQR(P*P+Q*Q)
900 LET C=SQR((V1+Q)/(2*V1))
910 LET S=P/(2*V1*C)
920 FOR I=1 TO M
930 LET V1=A(I,J)
940 LET A(I,J)=V1*C+A(I,K)*S
950 LET A(I,K)=-V1*S+A(I,K)*C
960 NEXT I
970 FOR I=1 TO N
980 LET V1=T(I,J)
990 LET T(I,J)=V1*C+T(I,K)*S
1000 LET T(I,K)=-V1*S+T(I,K)*C
1010 NEXT I
1020 LET N9=N2
1030 GOTO 1060
1040 LET N9=N9-1
1050 IF N9=0 THEN GOTO 1150
1051 REM ?? GOTO was EXIT for NS BASIC
1060 NEXT K
1070 NEXT J
1080 PRINT "SWEEP",I9,
1090 IF O1>0 THEN PRINT #O1,"SWEEP ",I9," ",
1100 IF 6*INT(I9/6)<>I9 THEN GOTO 680
1110 IF O1>0 THEN PRINT #O1
1120 IF I9>=30 THEN GOTO 1150
1130 PRINT 
1140 GOTO 680
1150 PRINT 
1160 IF O1>0 THEN PRINT #O1
1170 PRINT "CONVERGENCE AT SWEEP ",I9
1180 IF O1>0 THEN PRINT #O1,"CONVERGENCE AT SWEEP ",I9
1190 FOR J=1 TO N
1200 LET Q=0
1210 FOR I=1 TO M
1220 LET Q=Q+A(I,J)^2
1230 NEXT I
1240 LET Q=SQR(Q)
1250 IF Q=0 THEN GOTO 1290
1260 FOR I=1 TO M
1270 LET A(I,J)=A(I,J)/Q
1280 NEXT I
1290 LET Z(J)=Q
1300 NEXT J
1310 PRINT 
1320 PRINT "SINGULAR VALUES"
1340 FOR J=1 TO N
1350 PRINT Z(J),
1370 IF 5*INT(J/5)<>J THEN GOTO 1400
1380 PRINT 
1400 NEXT J
1410 PRINT 
1430 PRINT "VARIABLE # OF REGRESSAND",
1440 INPUT M2
1450 IF M2<=0 THEN GOTO 350
1470 LET S1=0
1480 FOR I=1 TO M
1490 LET S1=S1+(Y(I,M2)-E3*Y(M+1,M2))^2
1500 NEXT I
1510 FOR J=1 TO N
1520 LET S=0
1530 FOR I=1 TO M
1540 LET S=S+A(I,J)*Y(I,M2)
1550 NEXT I
1560 LET G(J)=S
1570 NEXT J
1580 PRINT "ENTER TOLERANCE FOR ZERO",
1590 INPUT Q
1600 IF Q<0 THEN GOTO 1410
1610 PRINT "SINGULAR VALUES <=",Q," ARE TAKEN AS 0"
1630 LET R=0
1640 FOR I=1 TO N
1650 LET V1=0
1660 LET S=0
1670 LET P=0
1680 FOR K=1 TO N
1690 LET C=0
1700 IF Z(K)<=Q THEN GOTO 1730
1710 LET C=1/Z(K)
1720 LET V1=V1+1
1730 LET S=S+C*T(I,K)*G(K)
1740 LET P=P+(C*T(I,K))^2
1750 NEXT K
1760 LET U(I)=P
1770 LET X(I)=S
1780 LET R=R+S*S
1790 NEXT I
1800 LET X(N)=X(N)*E3
1810 PRINT 
1820 PRINT "RESIDUALS"
1840 LET C=0
1850 LET S2=0
1860 FOR I=1 TO M
1870 LET S=Y(I,M2)-X(N)
1880 FOR K=1 TO N-1
1890 LET S=S-Y(I,W(K))*X(K)
1900 NEXT K
1910 PRINT S,
1930 IF 5*INT(I/5)<>I THEN GOTO 1960
1940 PRINT 
1960 LET C=C+S*S
1970 IF I=1 THEN GOTO 1990
1980 LET S2=S2+(S-S3)^2
1990 LET S3=S
2000 NEXT I
2010 PRINT 
2020 LET P=0
2040 IF M<=V1 THEN GOTO 2060
2050 LET P=C/(M-V1)
2060 PRINT M-V1," DEGREES OF FREEDOM"
2080 REM PRINT
2090 PRINT "SOLUTION VECTOR - CONSTANT LAST"
2110 FOR I=1 TO N
2120 LET V1=SQR(P*U(I))
2130 PRINT "X(",W(I),")=",X(I),"  STD.ERR.=",V1,
2140 IF O1>0 THEN PRINT #O1,"X(",W(I),")=",X(I),"  STD.ERR.=",V1,
2150 IF V1<=0 THEN GOTO 2180
2160 PRINT "  T=",ABS(X(I)/V1),
2170 IF O1>0 THEN PRINT #O1,"  T=",ABS(X(I)/V1),
2180 PRINT 
2190 IF O1>0 THEN PRINT #O1
2200 NEXT I
2210 PRINT "SUM OF SQUARES",C," SIGMA^2",P
2220 IF O1>0 THEN PRINT #O1,"SUM OF SQUARES",C," SIGMA^2",P
2230 PRINT "NORM OF SOLUTION",SQRT(R)
2240 IF O1>0 THEN PRINT #O1,"NORM OF SOLUTION",SQRT(R)
2250 PRINT "R SQUARED=",1-C/S1," DURBIN-WATSON STAT.=",S2/C
2260 IF O1>0 THEN PRINT #O1,"R SQUARED=",1-C/S1," DURBIN-WATSON STAT.=",S2/C
2270 PRINT 
2280 IF O1>0 THEN PRINT #O1
2290 GOTO 1580
2300 REM GET SERIES FROM FILE
2310 PRINT "FILENAME OR 'KEYBOARD' OR 'K'",
2320 INPUT G$
2330 IF LEN(G$)>0 THEN LET F$=G$ 
2331 REM DEFAULTS TO LAST SETTING
2340 PRINT "DATA FROM FILE :",F$
2350 IF F$="KEYBOARD" THEN 2420
2360 IF F$<>"K" THEN 2460
2370 PRINT
2380 PRINT "ENTER SERIES"
2390 FOR I=1 TO M
2400 INPUT1 Y(I,J)
2410 IF 5*INT(I/5)=I THEN PRINT
2420 NEXT I
2430 PRINT
2440 IF O1>0 THEN GOSUB 2860
2450 RETURN
2460 IF FILE(F$)=3 THEN 2490
2470 PRINT "FILE NOT FOUND OR OF WRONG TYPE"
2480 GOTO 2310
2490 OPEN #1,F$
2500 PRINT "SERIES NAME OR #",
2510 INPUT X$
2520 IF X$(1,1)="#" THEN 2770
2530 IF TYP(1)=0 THEN 2740
2540 IF TYP(1)=1 THEN 2570
2550 READ #1,C
2560 GOTO 2530
2570 READ #1,Y$
2580 IF X$<>Y$ THEN 2530
2590 I=0
2600 PRINT "SERIES:",Y$
2610 IF O1>0 THEN PRINT #O1,"SERIES:",Y$
2620 IF TYP(1)<>2 THEN 2690
2630 IF I=M THEN 2690
2640 I=I+1
2650 READ#1,Y(I,J)
2660 PRINT Y(I,J),
2670 IF 5*INT(I/5)=I THEN PRINT
2680 GOTO 2620
2690 PRINT
2700 PRINT "END OF SERIES  ",I," DATA POINTS"
2710 IF O1>0 THEN GOSUB 2860
2720 CLOSE #1
2730 RETURN
2740 PRINT "END OF FILE"
2750 CLOSE #1
2760 GOTO 2310
2770 X$=X$(2)
2780 P1=VAL(X$)
2790 J=0
2800 IF TYP(1)=0 THEN 2740
2810 IF TYP(1)=1 THEN 2840
2820 READ#1,C
2830 GOTO 2800
2840 J=J+1
2850 READ#1,Y$
2860 FOR I=1 TO M
2870 PRINT #O1,Y(I,J),
2880 IF 5*INT(I/5)=I THEN PRINT #O1
2890 NEXT I
2900 PRINT #O1
2910 RETURN
