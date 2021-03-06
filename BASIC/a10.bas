5 DIM A(10, 20),X(10),Y(10)
10 PRINT "GII JULY 25 77 ALG 10"
20 PRINT "GAUSS ELIMINATION FOR INVERSE ITERATION"
30 PRINT "ORDER=",
40 READ N
50 PRINT N
55 IF N <= 0 THEN QUIT : REM BWBASIC VARIANT
60 GOSUB 1500: REM BUILD OR INPUT MATRIX
70 GOSUB 2000: REM PUT METRIC IN RIGHT HALF OF A
75 GOSUB 1000: REM INITIAL GUESS TO VECTOR
80 LET K9=0 : REM SHIFT OF 0 FOR THIS EXAMPLE
90 PRINT "SHIFT=",K9
95 LET E9=K9
100 REM PRINT 
105 LET T2=N: REM FACTOR FOR CONVERGENCE TEST
110 LET T1=0: REM STEP 1
120 FOR I=1 TO N
130 LET Q=0
140 FOR J=1 TO N
150 LET A(I,J)=A(I,J)-K9*A(I,J+N)
160 LET S=S+ABS(A(I,J))
170 NEXT J
180 IF T1>=S THEN GOTO 200
190 LET T1=S
200 NEXT I
205 LET T1=T1*1.0E-7: REM NS 8 DIGIT BASIC
210 FOR I=1 TO N-1: REM STEP 2
218 LET S=ABS(A(I,I)): REM STEP 3
226 LET K=I
234 FOR J=I+1 TO N
242 IF ABS(A(J,I))<=S THEN GOTO 266
250 LET S=ABS(A(J,I))
258 LET K=J
266 NEXT J
274 IF S>0 THEN GOTO 298: REM STEP 4
282 LET A(I,I)=T1
290 GOTO 394
298 IF K=I THEN GOTO 346: REM STEP 5
306 FOR J=I TO 2*N: REM STEP 6
314 LET S=A(I,J)
322 LET A(I,J)=A(K,J)
330 LET A(K,J)=S
338 NEXT J
346 FOR J=I+1 TO N: REM STEP 7
354 LET S=A(J,I)/A(I,I)
362 FOR K=I TO 2*N
370 LET A(J,K)=A(J,K)-S*A(I,K)
378 NEXT K
386 NEXT J
394 NEXT I: REM STEP 8
402 IF ABS(A(N,N))>0 THEN GOTO 420: REM STEP 9
410 LET A(N,N)=T1
420 LET I9=0: REM STEP 10
430 LET I9=I9+1: REM STEP 11
440 LET M=N
445 LET S=X(N)
450 LET X(N)=Y(N)
455 LET Y(N)=S/A(N,N)
460 LET P=ABS(Y(N))
470 FOR I=(N-1) TO 1 STEP -1: REM STEP 12
480 LET S=X(I)
485 LET X(I)=Y(I)
490 FOR J=I+1 TO N
500 LET S=S-A(I,J)*Y(J)
510 NEXT J
520 LET Y(I)=S/A(I,I)
530 IF ABS(Y(I))<=P THEN GOTO 560
540 LET M=I
550 LET P=ABS(Y(I))
560 NEXT I
570 LET E8=K9+X(M)/Y(M): REM STEP 13
580 REM PRINT "APPROX EV=",E8
600 LET P=Y(M): REM STEP 14
610 LET M=0
620 FOR I=1 TO N
630 LET Y(I)=Y(I)/P
635 IF T2+Y(I)<>T2+X(I) THEN GOTO 640
636 LET M=M+1
640 NEXT I
644 IF M=N THEN GOTO 730: REM STEP 15 -- CONVERGENCE TEST
645 IF I9>100 THEN GOTO 730: REM LIMIT SET AT 100
650 FOR I=1 TO N: REM STEP 16
660 LET S=0
670 FOR J=1 TO N
680 LET S=S+A(I,J+N)*Y(J)
690 NEXT J
700 LET X(I)=S
710 NEXT I
720 GOTO 430: REM STEP 17
725 REM STEP 18 -- END AND RESIDUALS
730 PRINT "CONVERGED TO EV=",E8," IN ",I9," ITNS"
735 PRINT M," EQUAL CPNTS IN VECTOR BETWEEN ITERATIONS"
740 GOSUB 1500: REM GET MATRIX AGAIN
750 GOSUB 2000: REM GET METRIC AGAIN
755 LET S=0: REM COMPUTE VECTOR INNER PRODUCT
760 FOR I=1 TO N
770 FOR J=1 TO N
780 LET S=S+Y(I)*A(I,J+N)*Y(J)
790 NEXT J
800 NEXT I
810 LET S=1/SQR(S)
815 PRINT "VECTOR"
820 FOR I=1 TO N
830 LET Y(I)=Y(I)*S: REM VECTOR NORMALIZATION
840 PRINT Y(I);
845 IF I=5*INT(I/5) THEN PRINT
850 NEXT I
855 PRINT 
860 PRINT "RESIDUALS"
870 FOR I=1 TO N
880 LET S=0
890 FOR J=1 TO N: REM MATRIX * VECTOR - VALUE * METRIC * VECTOR
900 LET S=S+(A(I,J)-E8*A(I,J+N))*Y(J)
910 NEXT J
920 PRINT S;
925 IF 5*INT(I/5)=I THEN PRINT
930 NEXT I
940 PRINT 
950 GOTO 40 : REM NEXT TRY
960 DATA 5, 10, -1
1000 REM          INITIAL X
1010 FOR I=1 TO N
1020 LET X(I)=1: REM MAY BE A POOR CHOICE
1030 NEXT I
1040 RETURN 
1500 REM A IN FOR FRANK MATRIX
1505 PRINT "FRANK MATRIX"
1510 FOR I=1 TO N
1520 FOR J=1 TO I
1530 A(I,J)=J
1540 A(J,I)=J
1550 NEXT J
1560 NEXT I
1570 RETURN
2000 REM UNIT B IN RIGHT HALF OF MATRIX
2010 FOR I=1 TO N
2020 FOR J=1 TO N
2030 A(I,J+N)=0
2040 NEXT J
2050 A(I,I+N)=1
2060 NEXT I
2070 RETURN
