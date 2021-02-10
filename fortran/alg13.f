      SUBROUTINE A13ESV(N,A,NA,EPS,H,ISWP,IPR,Z)
C  ALGORITHM 13 EIGENPROBLEM OF A REAL SYMMETRIC MATRIX VIA SVD
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  N     =  ORDER OF PROBLEM
C  A     =  ARRAY CONTAINING MATRIX FOR WHICH EIGENVALUES ARE TO BE
C           COMPUTED.  RETURNS EIGENVECTORS AS COLUMNS
C  NA    =  FIRST DIMENSION OF A
C  EPS   =  MACHINE PRECISION
C  H     =  A NUMBER LARGER THAN ANY POSSIBLE EIGENVALUE. CHANGED
C           DURING EXECUTION. DO NOT ENTER AS A CONSTANT
C  ISWP  =  LIMIT ON SWEEPS (INPUT). SWEEPS USED (OUTPUT).
C  IPR   =  PRINT CHANNEL. IPR.GT.0 FOR PRINTING.
C  Z     =  EIGENVALUES (OUTPUT)
C  STEP 0
      INTEGER N,NA,ISWP,IPR,LISWP,I,J,COUNT,N1,J1
      REAL A(NA,N),EPS,H,V,Z(N),P,Q,R,S,C
      LISWP=ISWP
      ISWP=0
      N1=N-1
C  STEP 1
      DO 5 I=1,N
        V=A(I,I)
        DO 3 J=1,N
          IF(J.EQ.I)GOTO 3
          V=V-ABS(A(I,J))
   3    CONTINUE
        IF(V.LT.H)H=V
   5  CONTINUE
      IF(H.LE.EPS)GOTO 6
      H=0.0
      GOTO 30
   6  H=H-SQRT(EPS)
C  STEP 2
      DO 15 I=1,N
        A(I,I)=A(I,I)-H
  15  CONTINUE
C  STEP 3
  30  COUNT=0
C  CHECK FOR ORDER 1 PROBLEMS AND SKIP WORK
      IF(N.EQ.1)GOTO 160
      ISWP=ISWP+1
      IF(ISWP.GT.LISWP)GOTO 160
C  STEP 4
      DO 140 J=1,N1
C  STEP 5
        J1=J+1
        DO 130 K=J1,N
C  STEP 6
          P=0.0
          Q=0.0
          R=0.0
          DO 65 I=1,N
            P=P+A(I,J)*A(I,K)
            Q=Q+A(I,J)**2
            R=R+A(I,K)**2
  65      CONTINUE
C  STEP 7
          IF(1.0.LT.1.0+ABS(P/SQRT(Q*R)))GOTO 80
          IF(Q.LT.R)GOTO 80
          COUNT=COUNT+1
          GOTO 130
  80      Q=Q-R
C  STEP 8
          V=SQRT(4.0*P*P+Q*Q)
          IF(V.EQ.0.0)GOTO 130
C  STEP 9
          IF(Q.LT.0.0)GOTO 110
C  STEP 10
          C=SQRT((V+Q)/(2.0*V))
          S=P/(V*C)
          GOTO 120
C  STEP 11
 110      S=SQRT((V-Q)/(2.0*V))
          IF(P.LT.0.0)S=-S
          C=P/(V*S)
C  STEP 12
 120      DO 125 I=1,N
            V=A(I,J)
            A(I,J)=V*C+A(I,K)*S
            A(I,K)=-V*S+A(I,K)*C
 125      CONTINUE
C  STEP 13
 130    CONTINUE
C  STEP 14
 140  CONTINUE
C  STEP 15
      IF(IPR.GT.0)WRITE(IPR,970)ISWP,COUNT
 970  FORMAT( 9H AT SWEEP,I4,2X,I4,18H ROTATIONS SKIPPED)
      IF(COUNT.LT.N*(N-1)/2)GOTO 30
C  STEP 16
 160  DO 168 J=1,N
        S=0.0
        DO 162 I=1,N
          S=S+A(I,J)**2
 162    CONTINUE
        S=SQRT(S)
        DO 164 I=1,N
          A(I,J)=A(I,J)/S
 164    CONTINUE
        R=S+H
        Z(J)=R
 168  CONTINUE
C  STEP 17
 170  RETURN
      END

