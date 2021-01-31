      SUBROUTINE A3GR(M,N,A,NDIM,Q,EPS,SAVEQ)
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
      INTEGER N,M,NA,MN,I,J,K,J1
      REAL A(NDIM,N),Q(NDIM,M),EPS,TOL,B,P,S,C
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
          C=A(J,J)
          S=A(K,J)
          B=ABS(C)
          IF(ABS(S).GT.B)B=ABS(S)
          IF(B.EQ.0.0)GOTO 90
          C=C/B
          S=S/B
          P=SQRT(C*C+S*S)
C  STEP 4
          S=S/P
C  STEP 5
          IF(ABS(S).LT.TOL)GOTO 90
C  STEP 6
          C=C/P
C  STEP 7
          DO 75 I=1,N
            P=A(J,I)
            A(J,I)=C*P+S*A(K,I)
            A(K,I)=-S*P+C*A(K,I)
  75      CONTINUE
C  STEP 8
          IF(.NOT.SAVEQ)GOTO 90
          DO 85 I=1,M
            P=Q(I,J)
            Q(I,J)=C*P+S*Q(I,K)
            Q(I,K)=-S*P+C*Q(I,K)
  85      CONTINUE
C  STEP 9
  90    CONTINUE
C  STEP 10
 100  CONTINUE
      RETURN
      END

