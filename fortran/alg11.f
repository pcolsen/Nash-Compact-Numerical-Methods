      SUBROUTINE A11VS(N,X,Y,VNORM)
C  ALGORITHM 11  -  COMPLEX VECTOR STANDARDISATION
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  STANDARDISES COMPEX VECTOR (N ELEMENTS)  X+SQRT(-1)*Y
C      TO  1.0+SQRT(-1)*0.0
C  VNORM = NORM OF VECTOR (LARGEST ELEMENT)
C  STEP 0
      INTEGER N,I,K
      REAL X(N),Y(N),VNORM,G,B,E,S
C  STEP 1
      G=0.0
C  STEP 2
      DO 60 I=1,N
C  STEP 3
        B=X(I)**2+Y(I)**2
C  STEP 4
        IF(B.LE.G)GOTO 60
C  STEP 5
        K=I
        G=B
C  STEP 6
  60  CONTINUE
C  SAVE NORM
      VNORM=G
C  SAFETY CHECK
      IF(G.EQ.0.0)RETURN
C  STEP 7
      E=X(K)/G
      S=-Y(K)/G
C  STEP 8
      DO 85 I=1,N
        G=X(I)*E-Y(I)*S
        Y(I)=Y(I)*E+X(I)*S
        X(I)=G
  85  CONTINUE
C  END
      RETURN
      END
