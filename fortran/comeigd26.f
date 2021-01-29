C     MASTER COMDRIVE
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER UPP
      INTEGER CLKTC1, CLKTC2
      DIMENSION AR(10,10),AI(10,10),ASR(10,10),ASI(10,10),WR(10),WI(10)
      DIMENSION ZR(10,10),ZI(10,10),SCALE(10),INT(10)
      NIN = 5
      NOUT = 6
    1 READ (NIN,901) N
      IF (N.LE.0) STOP
      WRITE(NOUT, 900)N
 900  FORMAT(' COMPLEX MATRIX OF ORDER ',I4,' REAL FRANK, IMAG MOLER')
      DO 110 I=1,N
         DO 105 J=1,N
            AR(I,J) = MIN(I,J)
            AI(I,J) = MIN(I,J) - 2.0
  105    CONTINUE
         AR(I,I)=I
         AI(I,I)=I
  110  CONTINUE
C   START
      RGRD=0.0
      ISS=0
      JSS=0
      WRITE (NOUT,903) N
      WRITE (NOUT,904)
      do 103 i=1,n
      WRITE (NOUT,905) (AR(I,J),J=1,N)
 103  continue
      WRITE (NOUT,906)
      do 104 i=1,n
      WRITE (NOUT,905) (AI(I,J),J=1,N)
 104  continue
      DO 25 I=1,N
      DO 20 J=1,N
      ASR(I,J)=AR(I,J)
      ASI(I,J)=AI(I,J)
 20   CONTINUE
 25   CONTINUE
      RADX=16.0
c      CALL TIMER(CLKTC1)
c   20210128 -- SUPPRESS BALANCING FOR NOW
C      CALL CBAL(10,N,RADX,AR,AI,LOW,LUP,SCALE)
      CALL COMEIG(N,10,AR,AI,ZR,ZI,WR,WI)
C      CALL CBABK2(10,N,LOW,LUP,SCALE,N,ZR,ZI)
c      CALL TIMER(CLKTC2)
C      CLKTC2=CLKTC2-CLKTC1
      DO 6 J=1,N
      IT=1
      BIG=ZR(1,J)**2+ZI(1,J)**2
      IF (N.LT.2) GO TO 4
      DO 3 I=2,N
      U=ZR(I,J)**2+ZI(I,J)**2
      IF (U.LE.BIG) GO TO 3
      BIG=U
      IT=I
    3 CONTINUE
    4 U=ZR(IT,J)/BIG
      V=-ZI(IT,J)/BIG
      DO 5 I=1,N
      BIG=ZR(I,J)*U-ZI(I,J)*V
      ZI(I,J)=ZI(I,J)*U+ZR(I,J)*V
      ZR(I,J)=BIG
    5 CONTINUE
    6 CONTINUE
      DO 9 J=1,N
      WRITE (NOUT,907) J,WR(J),WI(J)
      WRITE (NOUT,908) (ZR(I,J),ZI(I,J),I=1,N)
      WRITE (NOUT,909)
      AL=WR(J)
      GA=WI(J)
      DO 8 I=1,N
      U=0.0
      V=0.0
      DO 7 K=1,N
      U=U+ASR(I,K)*ZR(K,J)-ASI(I,K)*ZI(K,J)
      V=V+ASR(I,K)*ZI(K,J)+ASI(I,K)*ZR(K,J)
    7 CONTINUE
      U=U-AL*ZR(I,J)+GA*ZI(I,J)
      V=V-GA*ZR(I,J)-AL*ZI(I,J)
      TEM=DSQRT(U**2+V**2)
      IF (TEM.LE.RGRD) GO TO 8
      RGRD=TEM
      ISS=I
      JSS=J
      WRITE (NOUT,908) U,V
    8 CONTINUE
    9 CONTINUE
      WRITE (NOUT,910) RGRD,JSS,ISS
C      WRITE (NOUT,912) CLKTC2
      GO TO 1
  901 FORMAT (I4)
  902 FORMAT (8F10.5)
  903 FORMAT (' ORDER OF MATRIX',I4)
  904 FORMAT (' REAL PART OF MATRIX')
  905 FORMAT (' ',1P5D16.8)
  906 FORMAT (' IMAGINARY PART OF MATRIX')
  907 FORMAT (' EIGENVALUE ',I4,' = ',1PD16.8,' + I*',1PD16.8)
  908 FORMAT (' ',2D16.8)
  909 FORMAT (' RESIDUALS,REAL AND IMAGINARY')
  910 FORMAT (' MAXIMUM RESIDUAL MAGNITUDE=',1PD16.8,' SOLUTION',I4,
     #' ELEMENT',I4)
C  912 FORMAT ('  TIME FOR EIGENSOLUTION =',I9,' * 0.01 SECS')
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     include 'comeig.for'
      SUBROUTINE COMEIG(N,ND,A,Z,T,U,RR,RI)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL MARK
C     SOLVES COMPLEX EIGENPROBLEM FOR A+I*Z
C     VECTORS (RIGHT HAND) RETURNED IN T+I*U
C     EIGENVALUES IN DECREASING ORDER OF MAGNITUDE DOWN DIAGONALS OF A
C     AND Z ARE RESTORED IN ER AND EI ON OUTPUT
C     EPS IS MACHINE DEPENDENT TOLERANCE
C     ITERATION LIMIT IS 35
      DIMENSION A(ND,N),Z(ND,N),T(ND,N),U(ND,N),RR(N),RI(N)
      EQUIVALENCE (AKI,AIK,TIK),(AIM,AMI,TIM)
      EQUIVALENCE (ZKI,ZIK,UIK),(ZMI,ZIM,UIM)
      IF (N.LE.1) GO TO 24
C     SET TOLERANCES
      CALL ENVROD(IB,IT,IR)
      EPS=(1.0D0*IB)**(1-IT)
      MARK=.FALSE.
C     PUT IDENTITY MATRIX IN T AND ZERO IN U
      N1=N-1
      DO 2 I=1,N1
      T(I,I)=1.0D0
      U(I,I)=0.0D0
      I1=I+1
      DO 1 J=I1,N
      T(I,J)=0.0D0
      U(I,J)=0.0D0
      U(J,I)=0.0D0
      T(J,I)=0.0D0
    1 CONTINUE
    2 CONTINUE
      T(N,N)=1.0D0
      U(N,N)=0.0D0
C     SAFETY LOOP
      DO 23 IT=1,35
      IF (MARK) GO TO 24
C     CONVERGENCE CRITERIA
      TAU=0.0D0
      DO 4 K=1,N
      TEM=0.0D0
      DO 3 I=1,N
      IF (I.NE.K) TEM=DABS(A(I,K))+DABS(Z(I,K))+TEM
    3 CONTINUE
      TAU=TAU+TEM
      RR(K)=TEM+DABS(A(K,K))+DABS(Z(K,K))
    4 CONTINUE
      WRITE (NOUT,901) TAU,IT
C     INTERCHANGE COLUMNS AND ROWS
      DO 8 K=1,N1
      SMAX=RR(K)
      I=K
      K1=K+1
      DO 5 J=K1,N
      IF (SMAX.GE.RR(J)) GO TO 5
      SMAX=RR(J)
      I=J
    5 CONTINUE
      IF (I.EQ.K) GO TO 8
      RR(I)=RR(K)
      DO 6 J=1,N
      TEP=A(K,J)
      A(K,J)=A(I,J)
      A(I,J)=TEP
      TEP=Z(K,J)
      Z(K,J)=Z(I,J)
      Z(I,J)=TEP
    6 CONTINUE
      DO 7 J=1,N
      TEP=A(J,K)
      A(J,K)=A(J,I)
      A(J,I)=TEP
      TEP=Z(J,K)
      Z(J,K)=Z(J,I)
      Z(J,I)=TEP
      TEP=T(J,K)
      T(J,K)=T(J,I)
      T(J,I)=TEP
      TEP=U(J,K)
      U(J,K)=U(J,I)
      U(J,I)=TEP
    7 CONTINUE
    8 CONTINUE
      IF (TAU.LT.(100.0D0*EPS)) GO TO 24
C     BEGIN SWEEP
      MARK=.TRUE.
      DO 22 K=1,N1
      K1=K+1
      DO 21 M=K1,N
      HJ=0.0D0
      HR=0.0D0
      HI=0.0D0
      G=0.0D0
      DO 9 I=1,N
      IF (I.EQ.K.OR.I.EQ.M) GO TO 9
      HR=HR+A(K,I)*A(M,I)+Z(K,I)*Z(M,I)-A(I,K)*A(I,M)-Z(I,K)*Z(I,M)
      HI=HI+Z(K,I)*A(M,I)-A(K,I)*Z(M,I)-A(I,K)*Z(I,M)+Z(I,K)*A(I,M)
      TE=A(I,K)**2+Z(I,K)**2+A(M,I)**2+Z(M,I)**2
      TEE=A(I,M)**2+Z(I,M)**2+A(K,I)**2+Z(K,I)**2
      G=G+TE+TEE
      HJ=HJ-TE+TEE
    9 CONTINUE
      BR=A(K,M)+A(M,K)
      BI=Z(K,M)+Z(M,K)
      ER=A(K,M)-A(M,K)
      EI=Z(K,M)-Z(M,K)
      DR=A(K,K)-A(M,M)
      DI=Z(K,K)-Z(M,M)
      TE=BR**2+EI**2+DR**2
      TEE=BI**2+ER**2+DI**2
      IF (TE.LT.TEE) GO TO 10
      SISW=1.0D0
      C=BR
      S=EI
      D=DR
      DE=DI
      ROOT2=DSQRT(TE)
      GO TO 11
   10 SISW=-1.0D0
      C=BI
      S=-ER
      D=DI
      DE=DR
      ROOT2=DSQRT(TEE)
   11 ROOT1=DSQRT(S*S+C*C)
      SIG=DSIGN(1.0D0,D)
      SA=0.0D0
      CA=DSIGN(1.0D0,C)
      IF (ROOT1.GE.EPS) GO TO 14
      SX=0.0D0
      SA=0.0D0
      CX=1.0D0
      CA=1.0D0
      IF (SISW.GT.0.0D0) GO TO 12
      E=EI
      B=-BR
      GO TO 13
   12 E=ER
      B=BI
   13 SND=D**2+DE**2
      GO TO 16
   14 IF (DABS(S).LE.EPS) GO TO 15
      CA=C/ROOT1
      SA=S/ROOT1
   15 COT2X=D/ROOT1
      COTX=COT2X+(SIG*DSQRT(1.0D0+COT2X**2))
      SX=SIG/DSQRT(1.0D0+COTX**2)
      CX=SX*COTX
C     FIND ROTATED ELEMENTS
      ETA=(ER*BR+BI*EI)/ROOT1
      TSE=(BR*BI-ER*EI)/ROOT1
      TE=SIG*(-ROOT1*DE+TSE*D)/ROOT2
      TEE=(D*DE+ROOT1*TSE)/ROOT2
      SND=ROOT2**2+TEE**2
      TEE=HJ*CX*SX
      COS2A=CA**2-SA**2
      SIN2A=2.0D0*CA*SA
      TEM=HR*COS2A+HI*SIN2A
      TEP=HI*COS2A-HR*SIN2A
      HR=CX*CX*HR-SX*SX*TEM-CA*TEE
      HI=CX*CX*HI+SX*SX*TEP-SA*TEE
      B=SISW*TE*CA+ETA*SA
      E=CA*ETA-SISW*TE*SA
   16 S=HR-SIG*ROOT2*E
      C=HI-SIG*ROOT2*B
      ROOT=DSQRT(C*C+S*S)
      IF (ROOT.GE.EPS) GO TO 17
      CB=1.0D0
      CH=1.0D0
      SB=0.0D0
      SH=0.0D0
      GO TO 18
   17 CB=-C/ROOT
      SB=S/ROOT
      TEE=CB*B-E*SB
      SNC=TEE*TEE
      TANH=ROOT/(G+2.0D0*(SNC+SND))
      CH=1.0D0/DSQRT(1.0D0-TANH**2)
      SH=CH*TANH
C     PREPARE FOR TRANSFORMATION
   18 TEM=SX*SH*(SA*CB-SB*CA)
      C1R=CX*CH-TEM
      C2R=CX*CH+TEM
      C1I=-SX*SH*(CA*CB+SA*SB)
      C2I=C1I
      TEP=SX*CH*CA
      TEM=CX*SH*SB
      S1R=TEP-TEM
      S2R=-TEP-TEM
      TEP=SX*CH*SA
      TEM=CX*SH*CB
      S1I=TEP+TEM
      S2I=TEP-TEM
C     DECIDE WHETHER TO MAKE TRANSFORMATION
      TEM=DSQRT(S1R**2+S1I**2)
      TEP=DSQRT(S2R**2+S2I**2)
      IF (TEM.LE.EPS.AND.TEP.LE.EPS) GO TO 21
      MARK=.FALSE.
C     TRANSFORMATION ON LEFT
      DO 19 I=1,N
      AKI=A(K,I)
      AMI=A(M,I)
      ZKI=Z(K,I)
      ZMI=Z(M,I)
      A(K,I)=C1R*AKI-C1I*ZKI+S1R*AMI-S1I*ZMI
      Z(K,I)=C1R*ZKI+C1I*AKI+S1R*ZMI+S1I*AMI
      A(M,I)=S2R*AKI-S2I*ZKI+C2R*AMI-C2I*ZMI
      Z(M,I)=S2R*ZKI+S2I*AKI+C2R*ZMI+C2I*AMI
   19 CONTINUE
C     TRANSFORMATION ON RIGHT
      DO 20 I=1,N
      AKI=A(I,K)
      AMI=A(I,M)
      ZKI=Z(I,K)
      ZMI=Z(I,M)
      A(I,K)=C2R*AKI-C2I*ZKI-S2R*AMI+S2I*ZMI
      Z(I,K)=C2R*ZKI+C2I*AKI-S2R*ZMI-S2I*AMI
      A(I,M)=-S1R*AKI+S1I*ZKI+C1R*AMI-C1I*ZMI
      Z(I,M)=-S1R*ZKI-S1I*AKI+C1R*ZMI+C1I*AMI
      AKI=T(I,K)
      AMI=T(I,M)
      ZKI=U(I,K)
      ZMI=U(I,M)
      T(I,K)=C2R*AKI-C2I*ZKI-S2R*AMI+S2I*ZMI
      U(I,K)=C2R*ZKI+C2I*AKI-S2R*ZMI-S2I*AMI
      T(I,M)=-S1R*AKI+S1I*ZKI+C1R*AMI-C1I*ZMI
      U(I,M)=-S1R*ZKI-S1I*AKI+C1R*ZMI+C1I*AMI
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
   23 CONTINUE
   24 DO 25 I=1,N
      RR(I)=A(I,I)
      RI(I)=Z(I,I)
   25 CONTINUE
      RETURN
  901 FORMAT (' TAU=',D20.10,' AT ITERATION ',I4)
      END
      SUBROUTINE ENVROD(BETA,T,RND)
C     DOUBLE PRECISION MACHINE ENVIRONMENT
C     NO INPUT
C     PARAMETERS:
C        BETA  - MACHINE RADIX
C        T     - NUMBER OF RADIX DIGITS IN WORD (DOUBLE)
C        RND   - SET TO 1 IF MACHINE ROUNDS
C                SET TO 0 IF MACHINE CHOPS
C     ALL PARAMETERS ARE INTEGERS
C     MACHINE DOUBLE PRECISION I.E. SMALLEST POSITIVE NUMBER SUCH
C     THAT 1. + NUMBER .GT. 1. , IS DBLE(FLOAT(BETA))**(1-T)
      INTEGER BETA,T,RND
      DOUBLE PRECISION A,B,DBLE
      RND=1
      A=2.0D0
      B=2.0D0
    1 IF ((A+1.0D0)-A.NE.1.0D0) GO TO 2
      A=2.0D0*A
      GO TO 1
    2 IF (A+B.NE.A) GO TO 3
      B=2.0D0*B
      GO TO 2
    3 BETA=(A+B)-A
      IF (A+DBLE(FLOAT(BETA-1)).EQ.A) RND=0
      IF (A+DBLE(FLOAT(BETA-1)).EQ.A) RND=0
      T=0
      A=1.0D0
    4 T=T+1
      A=A*DBLE(FLOAT(BETA))
      IF ((A+1.0D0)-A.EQ.1.0D0) GO TO 4
      RETURN
      END
C      include 'cbal.for'
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CBAL (NM,N,RADIX,AR,AI,LOW,UPP,SCALE)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      REAL*8 AR(NM,N),AI(NM,N),SCALE(N)
      INTEGER UPP
      LOGICAL NOCONV
      COMPLEX*16 DCMPLX
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBALANCE, WHICH IS A COMPLEX VERSION OF  BALANC,
C     NUM. MATH. 13,293-304(1969) BY PARLETT AND REINSCH.
C
C     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT--
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        RADIX IS THE BASE OF THE MACHINE FLOATING POINT
C          REPRESENTATION.  FOR SYSTEM/360 THIS IS 16.0,
C
C        AR AND AI ARE ARRAYS CONTAINING THE REAL AND IMAGINARY
C          PARTS, RESPECTIVELY, OF THE ELEMENTS OF THE MATRIX
C          TO BE BALANCED.
C
C     ON OUTPUT--
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE ELEMENTS OF THE BALANCED MATRIX,
C
C        LOW, UPP ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)
C          ARE EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=UPP+1,...,N,
C
C        SCALE IS AN ARRAY WHICH CONTAINS INFORMATION DETERMINING
C           THE PERMUTATIONS USED AND SCALING FACTORS.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH UPP
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J)       J = LOW,...,UPP
C                 = P(J)         J = UPP+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO UPP+1,
C     THEN 1 TO LOW-1.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN CBALANCE APPEARS IN
C     CBAL  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     ARITHMETIC IS REAL EXCEPT FOR THE USE OF THE SUBROUTINES
C     CABS AND CMPLX.
C
C     TRANSLATED BY V. KLEMA, ARGONNE NATIONAL LABORATORY, AUG., 1969.
C     MODIFIED BY B. GARBOW, JAN., 1971.
C
C     ----------------------------------------------------------------
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     ********** IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE **********
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN **********
   80 L = L - 1
      IF (L .LT. 1) GO TO 280
  100 LP1 = L + 1
C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********
      DO 120 JJ = 1, L
         J = LP1 - JJ
         R = 0.0
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            R = R + CDABS(DCMPLX(AR(J,I),AI(J,I)))
  110    CONTINUE
C
         IF (R .NE. 0.0) GO TO 120
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT **********
  130 K = K + 1
C
  140 DO 170 J = K, L
         C = 0.0
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            C = C + CDABS(DCMPLX(AR(I,J),AI(I,J)))
  150    CONTINUE
C
         IF (C .NE. 0.0) GO TO 170
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********
      DO 180 I = K, L
      SCALE(I) = 1.0
  180 CONTINUE
C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0
         R = 0.0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + CDABS(DCMPLX(AR(J,I),AI(J,I)))
            R = R + CDABS(DCMPLX(AR(I,J),AI(I,J)))
  200    CONTINUE
C
         G = R / RADIX
         F = 1.0
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     ********** NOW BALANCE **********
  240    IF (F .EQ. 1.0) GO TO 270
         G = 1.0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
C
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
C
  270 CONTINUE
C
      IF(NOCONV) GO TO 190
  280 LOW = K
      UPP = L
      RETURN
C     ********** LAST CARD OF CBAL **********
      END
C      include 'cbabk2.for'
C
C     ------------------------------------------------------------------
C
      SUBROUTINE CBABK2(NM,N,LOW,UPP,SCALE,M,ZR,ZI)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      REAL*8 SCALE(N),ZR(NM,M),ZI(NM,M)
      INTEGER UPP
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBABK2, WHICH IS A COMPLEX VERSION OF BALBAK,
C     NUM. MATH. 13,293-304(1969) BY PARLETT AND REINSCH.
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX
C     GENERAL MATRIX BY BACK TRANSFORMING THOSE OF THE
C     CORRESPONDING BALANCED MATRIX PRODUCED BY  CBAL.
C
C     ON INPUT--
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        LOW AND UPP ARE INTEGERS PRODUCED BY THE BALANCING
C          SUBROUTINE  CBAL,
C
C        SCALE IS AN ARRAY PRODUCED BY  CBAL  CONTAINING
C          INFORMATION ABOUT THE PERMUTATION AND SCALING
C          TRANSFORMATIONS USED IN BALANCING,
C
C        M IS AN INTEGER GIVING THE NUMBER OF COLUMNS OF
C          Z = (ZR,ZI) TO BE BACK TRANSFORMED,
C
C        ZR AND ZI ARE ARRAYS, THE FIRST M COLUMNS OF WHICH
C          CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C          OF THE EIGENVECTORS TO BE BACK TRANSFORMED.
C
C     ON OUTPUT--
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE BACK TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     TRANSLATED BY V. KLEMA, ARGONNE NATIONAL LABORATORY, AUG., 1969.
C     MODIFIED BY B. GARBOW, JAN., 1971.
C
C     ------------------------------------------------------------------
C
      IF (UPP .LT. LOW) GO TO 120
C
      DO 110 I = LOW, UPP
         S = SCALE(I)
C     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0/SCALE(I) **********
         DO 100 J = 1, M
            ZR(I,J) = ZR(I,J) * S
            ZI(I,J) = ZI(I,J) * S
  100    CONTINUE
C
  110 CONTINUE
C     ********** FOR I=LOW-1 STEP -1 UNTIL 1,
C                UPP+1 STEP 1 UNTIL N DO -- **********
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. UPP) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
      RETURN
C     ********** LAST CARD OF CBABK2 **********
      END
