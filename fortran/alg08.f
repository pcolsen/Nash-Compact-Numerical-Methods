      SUBROUTINE A8CS(A,N2,X,N)
C  ALGORITHM 8
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  CHOLESKI BACK-SOLUTION - ALGORITHM 8
C  STEP 0
      INTEGER N2,N,Q,I,I1,J,II,QJ
      REAL A(N2),X(N)
C  STEP 1
C  SAFETY CHECK ON N2
      IF(N2.NE.N*(N+1)/2)STOP
      IF(A(1).EQ.0.0)X(1)=0.0
      IF(A(1).GT.0.0)X(1)=X(1)/A(1)
C  STEP 2
      IF(N.EQ.1)GOTO 50
C  STEP 3
      Q=1
C  STEP 4
      DO 40 I=2,N
C  STEP 5
        I1=I-1
        DO 10 J=1,I1
          Q=Q+1
          X(I)=X(I)-A(Q)*X(J)
  10    CONTINUE
C  STEP 6
        Q=Q+1
C  STEP 7
        IF(A(Q).EQ.0.0)X(I)=0.0
        IF(A(Q).GT.0.0)X(I)=X(I)/A(Q)
C  STEP 8
  40  CONTINUE
C  STEP 9
  50  IF(A(N2).EQ.0.0)X(N)=0.0
      IF(A(N2).GT.0.0)X(N)=X(N)/A(N2)
C  STEP 10
      IF(N.EQ.1)GOTO 100
C  STEP 11
      DO 80 II=2,N
        I=N+2-II
C  STEP 12
        Q=I*(I-1)/2
C  STEP 13
        I1=I-1
        DO 60 J=1,I1
          QJ=Q+J
          X(J)=X(J)-X(I)*A(QJ)
  60    CONTINUE
C  STEP 14
        IF(A(Q).EQ.0.0)X(I1)=0.0
        IF(A(Q).GT.0.0)X(I1)=X(I1)/A(Q)
C  STEP 15
  80  CONTINUE
C  STEP 16
 100  RETURN
      END

