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
