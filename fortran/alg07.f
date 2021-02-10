      SUBROUTINE A7CH(A,N2,N,INDEF)
C  ALGORITHM 7
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  CHOLESKI DECOMPOSITION OF REAL-SYMMETRIC
      LOGICAL INDEF
      INTEGER N2,N,I,J,Q,M,K,J1,MK,QK
      REAL A(N2),S
      INDEF=.FALSE.
C  STEP 1
      DO 100 J=1,N
C  STEP 2
        Q=J*(J+1)/2
C  STEP 3
        IF(J.EQ.1)GOTO 50
C  STEP 4
        DO 40 I=J,N
          M=I*(I-1)/2+J
          S=A(M)
          J1=J-1
          DO 20 K=1,J1
            MK=M-K
            QK=Q-K
            S=S-A(MK)*A(QK)
  20      CONTINUE
          A(M)=S
  40    CONTINUE
C  STEP 5
  50    IF(A(Q).GT.0.0)GOTO 60
C  SET FLAG IN THIS CASE
        INDEF=.TRUE.
C  STEP 6
        A(Q)=0.0
C  ASSUMES MATRIX NON-NEGATIVE DEFINITE
C  STEP 7
  60    S=SQRT(A(Q))
C  STEP 8
        DO 80 I=J,N
          M=I*(I-1)/2+J
          IF(S.EQ.0.0)A(M)=0.0
          IF(S.GT.0.0)A(M)=A(M)/S
  80    CONTINUE
C  STEP 9
 100  CONTINUE
      RETURN
      END

