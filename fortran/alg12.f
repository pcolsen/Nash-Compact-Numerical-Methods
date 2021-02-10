      SUBROUTINE A12CVR(N,X,Y,A,NA,Z,NZ,E,F,U,V)
C  ALGORITHM 12 RESIDUALS OF A COMPLEX EIGENSOLUTION
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  N    =  ORDER OF PROBLEM
C  U + I*V =  ( A + I*Z - E - I*F)*(X + I*Y)  WHERE I=SQRT(-1)
C  NA,NZ = FIRST DIMENSIONS OF A & Z RESPECTIVELY
C  STEP 0
      INTEGER N,NA,NZ,J,K
      REAL A(NA,N),Z(NZ,N),X(N),Y(N),E,F,U(N),V(N),S,G
C  STEP 1
      DO 50 J=1,N
C  STEP 2
        S=-E*X(J)+F*Y(J)
        G=-F*X(J)-E*Y(J)
C  STEP 3
        DO 35 K=1,N
          S=S+A(J,K)*X(K)-Z(J,K)*Y(K)
          G=G+A(J,K)*Y(K)+Z(J,K)*X(K)
  35    CONTINUE
C  STEP 4  NOTE SAVE IN U & V
        U(J)=S
        V(J)=G
C  STEP 5
  50  CONTINUE
      RETURN
      END
