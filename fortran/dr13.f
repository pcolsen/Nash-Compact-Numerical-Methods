C&&& A13
C  TEST ALG. 13    JULY 1978
C  J.C. NASH   JULY 1978, APRIL 1989
      REAL H,EPS
      INTEGER N,ND,I,NOUT,NIN
      REAL A(10,10),B(10,10),AT(10,10),Z(10),V(10,10),RMAX,VMAX
      EXTERNAL FRANKM,UNITM
      ND=10
C  I/O CHANNELS
      NIN=5
      NOUT=6
   1  READ(NIN,900)N
 900  FORMAT(I4)
      WRITE(NOUT,901)N
 901  FORMAT(' ORDER N=',I4)
      IF(N.LE.0)STOP
      CALL FRANKM(N,N,V,ND)
      ISWP=30
C   IBM  SHORT PRECISION
      EPS=16.0**(-5)
C  IBM VALUE FOR BIG NO.
C&&&       H=R1MACH(2)
      H = 1.0E+35
      CALL A13ESV(N,V,ND,EPS,H,ISWP,NOUT,Z)
      WRITE(NOUT,903)ISWP
 903  FORMAT(' CONVERGED IN ',I4,' SWEEPS')
      CALL EVT(N,V,ND,Z,FRANKM,UNITM,AT,ND,B,ND,NOUT,RMAX,VMAX)
      GOTO 1
      END
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
      SUBROUTINE UNITM(M,N,A,NA)
C  PUTS UNIT MATRIX M BY N IN A
C  J.C. NASH   JULY 1978, APRIL 1989
      INTEGER M,N,NA,I,J
      REAL A(NA,N)
      DO 10 I=1,M
        DO 5 J=1,N
          A(I,J)=0.0
          IF(I.EQ.J)A(I,I)=1.0
   5    CONTINUE
  10  CONTINUE
      RETURN
      END
      SUBROUTINE EVT(N,V,NV,Z,AIN,BIN,A,NA,B,NB,NOUT,RMAX,VMAX)
C  J.C. NASH   JULY 1978, APRIL 1989
C  COMPUTES RESIDUALS AND INNER PRODUCTS
C    R  =  (A - Z(J)*B)*V(.,J)
C  AIN AND BIN ARE NAMES OF MATRIX CALCULATING ROUTINES FOR A AND B
C  WHOSE FIRST DIMENSIONS ARE NA AND NB RESP.
C  RMAX AND VMAX ARE MAX ABS RESIDUAL AND INNER PRODUCT RESP.
C
      INTEGER N,NV,NA,NB,NOUT,I,J,K,RPOSI,RPOSJ,VPOSI,VPOSJ,I1,N1
      REAL V(NV,N),A(NA,N),B(NB,N),Z(N),RMAX,VMAX
      DOUBLE PRECISION ACC,TACC,DABS,DBLE
      CALL AIN(N,N,A,NA)
      CALL BIN(N,N,B,NB)
      N1=N-1
      TACC=0.0
      RPOSI=1
      RPOSJ=1
      DO 20 I=1,N
        DO 15 J=1,N
          ACC=0.0
          DO 10 K=1,N
            ACC=ACC+DBLE(V(K,J))*(A(I,K)-Z(J)*B(I,K))
  10      CONTINUE
          IF(DABS(ACC).LE.TACC)GOTO 15
          TACC=DABS(ACC)
          RPOSI=I
          RPOSJ=J
  15    CONTINUE
  20  CONTINUE
      RMAX=TACC
      IF(NOUT.GT.0)WRITE(NOUT,951)RMAX,RPOSI,RPOSJ
 951  FORMAT(' MAX. ABS. RESIDUAL=',1PE16.8,'  POSN',2I4)
      VPOSI=0
      VPOSJ=0
      TACC=0.0
      IF(N.EQ.1)GOTO 45
      DO 40 I=1,N1
        I1=I+1
        DO 35 J=I1,N
          ACC=0.0
          DO 30 K=1,N
            ACC=ACC+DBLE(V(K,I))*V(K,J)
  30      CONTINUE
          IF(DABS(ACC).LE.TACC)GOTO 35
          TACC=DABS(ACC)
          VPOSI=I
          VPOSJ=J
  35    CONTINUE
  40  CONTINUE
      VMAX=TACC
      IF(NOUT.GT.0)WRITE(NOUT,952)VMAX,VPOSI,VPOSJ
 952  FORMAT(' MAX. ABS. INNER PRODUCT=',1PE16.8,'  POSN',2I4)
  45  IF(NOUT.LE.0)RETURN
      DO 50 J=1,N
        WRITE(NOUT,953)J,Z(J)
 953  FORMAT(' EIGENVALUE',I3,'=',1PE16.8)
        WRITE(NOUT,954)(V(K,J),K=1,N)
 954  FORMAT(1H ,5E16.8)
  50  CONTINUE
      RETURN
      END
      SUBROUTINE FRANKM(M,N,A,NA)
C  J.C. NASH   JULY 1978, APRIL 1989
      INTEGER M,N,NA,I,J
C  INPUTS FRANK MATRIX M BY N INTO A
      REAL A(NA,N)
      DO 20 I=1,M
        DO 10 J=1,N
          A(I,J)=AMIN0(I,J)
  10    CONTINUE
  20  CONTINUE
      RETURN
      END
