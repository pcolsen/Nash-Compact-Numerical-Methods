C&&& A7-8
C  TEST ALGORITHMS 7 & 8 -- A7CH +A8CS
C  J.C. NASH   JULY 1978, APRIL 1989
C  USES FRANK MATRIX
      LOGICAL INDEF
      INTEGER N,N2,I,J,IJ,NOUT
      REAL A(55),G(10)
C  PRINTER CHANNEL
      NOUT=6
C  MAIN LOOP
      DO 100 N=2,10,2
        N2=N*(N+1)/2
        WRITE(NOUT,950)N
 950    FORMAT('0ORDER=',I3,'  ORIGINAL FRANK MATRIX')
        DO 20 I=1,N
          DO 10 J=1,I
            IJ=I*(I-1)/2+J
            A(IJ)=J
  10      CONTINUE
  20    CONTINUE
        CALL SOUT(A,N2,N,NOUT)
        CALL A7CH(A,N2,N,INDEF)
      WRITE(NOUT,955)
 955  FORMAT('0DECOMPOSED MATRIX')
      CALL SOUT(A,N2,N,NOUT)
        IF(INDEF)WRITE(NOUT,953)
 953    FORMAT('0INDEFINITE MATRIX')
C  COMPUTE RHS
      DO 30 I=1,N
        G(I)=1.0
 30   CONTINUE
      WRITE(NOUT,954)(G(J),J=1,N)
 954  FORMAT(' RHS',1P5E16.8)
      CALL A8CS(A,N2,G,N)
      WRITE(NOUT,957)(G(J),J=1,N)
 957  FORMAT(' SOL',1P5E16.8)
 100  CONTINUE
      STOP
      END
      SUBROUTINE SOUT(A,N2,N,NOUT)
C  J.C. NASH   JULY 1978, APRIL 1989
      INTEGER N2,N,NOUT,I,J,IJ,JJ
      REAL A(N2)
C     PRINTS SYMMETRIC MATRIX STORED ROW-WISE AS A VECTOR
      DO 20 I=1,N
        WRITE(NOUT,951)I
 951    FORMAT(' ROW',I3)
        IJ=I*(I-1)/2+1
        JJ=IJ+I-1
        WRITE(NOUT,952)(A(J),J=IJ,JJ)
 952    FORMAT(1H ,1P5E16.8)
 20   CONTINUE
      RETURN
      END
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
