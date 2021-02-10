      SUBROUTINE A9GJ(A,N2,N,INDEF,X)
C  ALGORITHM 9
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  BAUER-REINSCH  GAUSS-JORDAN INVERSION OF A SYMMETRIC, POSITIVE
C  A=MATRIX - STORED AS A VECTOR -- ELEMENT I,J IN POSITION I*(I-1)/2+J
C  N2=LENGTH OF VECTOR A = N*(N+1)/2
C  N=ORDER OF MATRIX
C  INDEF=LOGICAL FLAG SET .TRUE. IF MATRIX NOT COMPUTATIONALLY
C     POSITIVE DEFINITE
C  X=WORKING VECTOR OF LENGTH AT LEAST N
C  DEFINITE MATRIX
C  STEP 0
      LOGICAL INDEF
      INTEGER N2,N,K,KK,Q,M,Q2,JI,JQ
      REAL A(N2),S,T,X(N)
C  STEP 1
      INDEF=.FALSE.
      DO 100 KK=1,N
        K=N+1-KK
C  STEP 2
        S=A(1)
C  STEP 3
        IF(S.LE.0.0)INDEF=.TRUE.
        IF(INDEF)RETURN
C  STEP 4
        M=1
C  STEP 5
        DO 60 I=2,N
C  STEP 6
          Q=M
          M=M+I
          T=A(Q+1)
          X(I)=-T/S
C  STEP 7
          Q2=Q+2
          IF(I.GT.K)X(I)=-X(I)
C  STEP 8
          DO 40 J=Q2,M
            JI=J-I
            JQ=J-Q
            A(JI)=A(J)+T*X(JQ)
  40      CONTINUE
C  STEP 9
  60    CONTINUE
C  STEP 10
        Q=Q-1
        A(M)=1/S
C  STEP 11
        DO 80 I=2,N
          JI=Q+I
          A(JI)=X(I)
  80    CONTINUE
C  STEP 12
 100  CONTINUE
      RETURN
      END

