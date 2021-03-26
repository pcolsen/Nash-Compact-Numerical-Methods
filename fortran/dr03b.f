C&&& A3
C  TEST ALGORITHM 3
C  J.C. NASH   JULY 1978, APRIL 1989
      LOGICAL SAVEQ
      CHARACTER QSAVE
      INTEGER M,N,NIN,NOUT
      DOUBLE PRECISION A(10,10),Q(10,10),EPS,S,W(10,10)
      NDIM=10
C  I/O CHANNELS
      NIN=5
      NOUT=6
   1  READ(NIN,900)M,N,QSAVE
 900  FORMAT(2I5,A1)
      WRITE(NOUT,950)M,N,QSAVE
 950  FORMAT('M=',I5,'  N=',I5,'  QSAVE=',A1)
      IF(M.EQ.0.OR.N.EQ.0)STOP
      SAVEQ=.FALSE.
      IF (QSAVE .EQ. "T") SAVEQ=.TRUE.
      CALL FRANKM(M,N,A,10)
      WRITE(NOUT,952)
 952  FORMAT('INITIAL MATRIX')
      CALL OUT(A,NDIM,M,N,NOUT)
      DO 10 I=1,M
        DO 5 J=1,N
C         COPY MATRIX TO WORKING ARRAY
          W(I,J)=A(I,J)
  5     CONTINUE
 10   CONTINUE
C  IBM MACHINE PRECISION
!      EPS=16.0**(-5)
      CALL A3GR(M,N,W,10,Q,SAVEQ)
      WRITE(NOUT,953)
 953  FORMAT('FULL DECOMPOSED MATRIX')
      CALL OUT(A,NDIM,M,N,NOUT)
      IF(SAVEQ)CALL A3DT(M,N,W,NDIM,Q,NOUT,A)
      GOTO 1
      END
C-------------------------------------
      SUBROUTINE A3DT(M,N,W,NDIM,Q,NOUT,A)
C  TESTS GIVENS' DECOMPOSITION
C  J.C. NASH   JULY 1978, APRIL 1989
      INTEGER M,N,NDIM,NOUT,I,J,K
      DOUBLE PRECISION A(NDIM,N),Q(NDIM,M),W(NDIM,N),S,T
      WRITE(NOUT,960)
 960  FORMAT(' Q MATRIX')
      CALL OUT(Q,NDIM,M,M,NOUT)
      WRITE(NOUT,961)
 961  FORMAT(' R MATRIX (STORED IN W')
      CALL OUT(W,NDIM,M,N,NOUT)
      IF(N.LT.M)GOTO 9
      S=1.0
      DO 5 I=1,M
        S=S*W(I,I)
   5  CONTINUE
      WRITE(NOUT,963)S
 963  FORMAT(' DETERMINANT=',1PE16.8)
   9  CONTINUE
      T=0.0
      DO 20 I=1,M
        DO 15 J=1,N
          S=0.0
          DO 10 K=1,M
            S=S+Q(I,K)*W(K,J)
  10      CONTINUE
          S=S-A(I,J)
          IF(ABS(S).GT.T)T=ABS(S)
  15    CONTINUE
  20  CONTINUE
      WRITE(NOUT,962)T
 962  FORMAT(' MAX. DEVN. OF RECONSTRUCTION FROM ORIGINAL=',E16.8)
      RETURN
      END
C-------------------------------------
      SUBROUTINE OUT(A,NDIM,N,NP,NOUT)
C  J.C. NASH   JULY 1978, APRIL 1989
      INTEGER NDIM,N,NOUT,I,J
      DOUBLE PRECISION A(NDIM,NP)
      DO 20 I=1,N
        WRITE(NOUT,951)I
 951    FORMAT(' ROW',I3)
        WRITE(NOUT,952)(A(I,J),J=1,NP)
 952    FORMAT(1H ,1P5E16.8)
  20  CONTINUE
      RETURN
      END
C-------------------------------------
      SUBROUTINE FRANKM(M,N,A,NA)
C  J.C. NASH   JULY 1978, APRIL 1989
      INTEGER M,N,NA,I,J
C  INPUTS FRANK MATRIX M BY N INTO A
      DOUBLE PRECISION A(NA,N)
      DO 20 I=1,M
        DO 10 J=1,N
          A(I,J)=AMIN0(I,J)
  10    CONTINUE
  20  CONTINUE
      RETURN
      END
C==========================================
      SUBROUTINE A3GR(M,N,A,NDIM,Q,SAVEQ)
C Amended variant to use ideas of @Bindel2002Givens
C  ALGORITHM 3  GIVENS' REDUCTION
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  M,N  =  ORDER OF MATRIX TO BE DECOMPOSED
C  A    =  ARRAY CONTAINING MATRIX TO BE DECOMPOSED
C  NDIM   =  FIRST DIMENSION OF MATRICES - NDIM.GE.M
C  Q    =  ARRAY CONTAINING ORTHOGONAL MATRIX OF ACCUMULATED ROTATIONS
C  EPS  =  MACHINE PRECISION  = SMALLEST NO.GT.0.0 S.T. 1.0+EPS.GT.1.0
C  SAVEQ=  LOGICAL FLAG SET .TRUE. IF Q TO BE FORMED
C  STEP 0
      LOGICAL SAVEQ
      INTEGER N,M,NA,MN,I,J,K,J1,COUNT
      DOUBLE PRECISION A(NDIM,N),Q(NDIM,M),TOL,B,P,S,C,F,G,F1,G1,RC
      DOUBLE PRECISION ONE, TWO, ZERO
      PARAMETER       ( ONE = 1.0D+0, ZERO = 0.0D+0, TWO = 2.0D+0 )
      DOUBLE PRECISION SAFMIN, EPS, SAFEPS, SAFMN2, SAFMX2, LSE
      DOUBLE PRECISION LL2, LL22, ASFMN2
      SAFMIN = 2.2250738585072014D-308 ! DLAMCH( 'S' )
      EPS = TWO**(-53) ! DLAMCH( 'E' )
      SAFEPS=SAFMIN/EPS
      LSE=LOG(SAFEPS)
      LL2=LSE/LOG(TWO)
      LL22=LL2/TWO
      SAFMN2 = TWO**INT(LL22)
      SAFMX2 = ONE / SAFMN2

      write(*,*)"SAFMN2=",SAFMN2,"  SAFMX2=",SAFMX2

!      ROTCHK = 0.05 ! used for inconsequential rotation test. Was 0.25
      MN=M
      IF(M.GT.N)MN=N
      IF(.NOT.SAVEQ)GOTO 9
      DO 5 I=1,M
        DO 4 J=1,M
          Q(I,J)=0.0
   4    CONTINUE
        Q(I,I)=1.0
   5  CONTINUE
   9  TOL=EPS*EPS
C  STEP 1
      IF(M.EQ.1)RETURN
      DO 100 J=1,MN
        J1=J+1
        IF(J1.GT.M)GOTO 100
C  STEP 2
        DO 90 K=J1,M
C  STEP 3
          G=A(K,J)
          F=A(J,J)
          IF( G.EQ.ZERO ) THEN
              C = ONE
              S = ZERO
              RC = F
          ELSE IF( F.EQ.ZERO ) THEN
              C = ZERO
              S = ONE
              RC = G
          ELSE
              F1 = F
              G1 = G
              SCALE = MAX( ABS( F1 ), ABS( G1 ) )
              IF( SCALE.GE.SAFMX2 ) THEN
              COUNT = 0
510         CONTINUE
              COUNT = COUNT + 1
              F1 = F1*SAFMN2
              G1 = G1*SAFMN2
              SCALE = MAX( ABS( F1 ), ABS( G1 ) )
              IF( SCALE.GE.SAFMX2 .AND. COUNT .LT. 20) GO TO 510
              RC = SQRT( F1**2+G1**2 )
              C = F1 / RC
              S = G1 / RC
              DO II = 1,COUNT
                RC = RC*SAFMX2
 520          END DO
            ELSE IF( SCALE.LE.SAFMN2 ) THEN
              COUNT = 0
 530          CONTINUE
              COUNT = COUNT + 1
              F1 = F1*SAFMX2
              G1 = G1*SAFMX2
              SCALE = MAX( ABS( F1 ), ABS( G1 ) )
              IF( SCALE.LE.SAFMN2 ) GO TO 530
              RC = SQRT( F1**2+G1**2 )
              C = F1 / RC
              S = G1 / RC
              DO II = 1, COUNT
                RC = RC*SAFMN2
              END DO
            ELSE
              RC = SQRT( F1**2+G1**2 )
              C = F1 / RC
              S = G1 / RC
            END IF
            IF( ABS( F ).GT.ABS( G ) .AND. C.LT.ZERO ) THEN
              C = -C
              S = -S
              RC = -RC
            END IF
          END IF
C  STEP 7
          DO 75 I=1,N
            P=A(J,I)
            A(J,I)=C*P+S*A(K,I)
            A(K,I)=-S*P+C*A(K,I)
  75      CONTINUE
C  STEP 8
          IF(.NOT.SAVEQ)GOTO 90
          DO I=1,M
            P=Q(I,J)
            Q(I,J)=C*P+S*Q(I,K)
            Q(I,K)=-S*P+C*Q(I,K)
          END DO
C  STEP 9
  90    CONTINUE
C  STEP 10
 100  CONTINUE
      RETURN
      END
