      SUBROUTINE A1SVD(M,N,A,NA,EPS,V,NV,Z,IPR)
C  ALGORITHM 1  SINGULAR VALUE DECOMPOSITION BY COLUMN ORTHOGONA-
C     LISATION VIA PLANE ROTATIONS
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  M BY N  MATRIX A  IS DECOMPOSED TO  U*Z*VT
C   A    =   ARRAY CONTAINING A (INPUT),  U (OUTPUT)
C   NA   =   FIRST DIMENSION OF A
C   EPS  =   MACHINE PRECISION
C   V    =   ARRAY IN WHICH ORTHOGAONAL MATRIX V IS ACCUMULATED
C   NV   =   FIRST DIMENSION OF V
C   Z    =   VECTOR OF SINGULAR VALUES
C  IPR   =  PRINT CHANNEL   IF IPR.GT.0  THEN PRINTING
C  STEP 0
      INTEGER M,N,J1,N1,COUNT
      REAL A(NA,N),V(NV,N),Z(N),EPS,TOL,P,Q,R,VV,C,S
C  UNDERFLOW AVOIDANCE STRATEGY
      REAL SMALL
      SMALL=1.0E-36
C  ABOVE IS VALUE FOR IBM
      TOL=N*N*EPS*EPS
      DO 6 I=1,N
        DO 4 J=1,N
          V(I,J)=0.0
   4    CONTINUE
      V(I,I)=1.0
   6  CONTINUE
      N1=N-1
C  STEP 1
  10  COUNT=N*(N-1)/2
C  STEP 2
      DO 140 J=1,N1
C  STEP 3
        J1=J+1
        DO 130 K=J1,N
C  STEP 4
          P=0.0
          Q=0.0
          R=0.0
C  STEP 5
          DO 55 I=1,M
            IF(ABS(A(I,J)).GT.SMALL.AND.ABS(A(I,K)).GT.SMALL)
     #         P=P+A(I,J)*A(I,K)
            IF(ABS(A(I,J)).GT.SMALL)Q=Q+A(I,J)**2
            IF(ABS(A(I,K)).GT.SMALL)R=R+A(I,K)**2
C           P=P+A(I,J)*A(I,K)
C           Q=Q+A(I,J)**2
C           R=R+A(I,K)**2
  55      CONTINUE
C  STEP 6
          IF(Q.GE.R)GOTO 70
          C=0.0
          S=1.0
          GOTO 90
C  STEP 7
  70      IF(R.LE.TOL)GOTO 120
          IF((P*P)/(Q*R).LT.TOL)GOTO 120
C  STEP 8
          Q=Q-R
          VV=SQRT(4.0*P**2+Q**2)
          C=SQRT((VV+Q)/(2.0*VV))
          S=P/(VV*C)
C  STEP 9
  90      DO 95 I=1,M
            R=A(I,J)
            A(I,J)=R*C+A(I,K)*S
            A(I,K)=-R*S+A(I,K)*C
  95      CONTINUE
C  STEP 10
          DO 105 I=1,N
            R=V(I,J)
            V(I,J)=R*C+V(I,K)*S
            V(I,K)=-R*S+V(I,K)*C
 105      CONTINUE
C  STEP 11
          GOTO 130
 120      COUNT=COUNT-1
C  STEP 13
 130    CONTINUE
C  STEP 14
 140  CONTINUE
C  STEP 15
      IF(IPR.GT.0)WRITE(IPR,964)COUNT
 964  FORMAT(1H ,I4,10H ROTATIONS)
      IF(COUNT.GT.0)GOTO 10
C  STEP 16
      DO 220 J=1,N
C  STEP 17
        Q=0.0
C  STEP 18
        DO 185 I=1,M
            Q=Q+A(I,J)**2
 185    CONTINUE
C  STEP 19
        Q=SQRT(Q)
        Z(J)=Q
        IF(IPR.GT.0)WRITE(IPR,965)J,Q
 965    FORMAT( 4H SV(,I3,2H)=,1PE16.8)
C  STEP 20
        IF(Q.LT.TOL)GOTO 220
C  STEP 21
        DO 215 I=1,M
          A(I,J)=A(I,J)/Q
 215    CONTINUE
C  STEP 22
 220  CONTINUE
      RETURN
      END

