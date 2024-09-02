C&&& A25
C  TEST ALG 25 USING GRID (5 POINT)
C  J.C. NASH   JULY 1978, APRIL 1989
      LOGICAL IFR
      INTEGER N,M,NOUT,NIN,KPR,LIMIT,I
      EXTERNAL APR,BPR
C     REAL EPS,P0,X(N),S(N),T(N),U(N),V(N),W(N),Y(N),RNORM
      COMMON /GSZ/ M,IFR,R(1600)
      REAL EPS,P0,RNORM,VNORM,RNV
      REAL S(1600),T(1600),U(1600),V(1600),W(1600),X(1600),Y(1600)
C  I/O CHANNELS
      NIN=5
      NOUT=6
   1  READ(NIN,900)N,LIMIT
 900  FORMAT(2I4)
      WRITE(NOUT,950)N,LIMIT
 950  FORMAT(' MOLER MATRIX ORDER',I5,'  LIMIT=',I4)
      IF(M.LE.0)STOP
      IFR=.FALSE.
C APPROX
      EPS=16.0**(-14) 
      KPR=LIMIT
      RNORM=1.0/SQRT(FLOAT(N))
      DO 10 I=1,N
        X(I)=RNORM
  10  CONTINUE
      CALL A25RQM(N,X,EPS,KPR,S,T,U,V,W,Y,P0,NOUT,APR,BPR)
      WRITE(NOUT,951)KPR,P0
 951  FORMAT(' RETURNED AFTER',I4,' PRODUCTS WITH EV=',1PE16.8)
      DO 20 I=1,N
        R(I)=-P0*X(I)
  20  CONTINUE
      CALL APR(N,X,V)
      RNORM=0.0
      VNORM=0.0
      DO 30 I=1,N
        RNORM=RNORM+(V(I)+R(I))**2
        VNORM=VNORM+X(I)**2
  30  CONTINUE
      RNORM=SQRT(RNORM/N)
      VNORM=SQRT(VNORM/N)
      RNV=RNORM/VNORM
      WRITE(NOUT,952)RNORM,VNORM,RNV
 952  FORMAT(' RESIDUAL NORM=',1PE16.8,' /',E16.8,'=',E16.8)
      GOTO 1
      END
      SUBROUTINE BPR(N,X,V)
C  J.C. NASH   JULY 1978, APRIL 1989
C  UNITM MATRIX * X  INTO V
      INTEGER N,I
      REAL X(N),V(N)
      DO 100 I=1,N
        V(I)=X(I)
 100  CONTINUE
      RETURN
      END
      SUBROUTINE APR(N,X,V)
      integer n, i, j
      double precision x(n), ax(n), sum
c     return ax = A * x for A = moler matrix
c     A[i,j]=min(i,j)-2 for i<>j, or i for i==j
      do 20 i=1,n
         sum=0.0
         do 10 j=1,n
            if (i.eq.j) then
               sum = sum+i*x(i)
            else
               sum = sum+(min(i,j)-2)*x(j)
            endif
 10      continue
         V(i)=sum
 20   continue
      return
      end
      RETURN
      END
      SUBROUTINE A25RQM(N,X,EPS,KPR,Y,Z,T,G,A,B,P0,IPR,APR,BPR)
C  ALGORITHM 25 RAYLEIGH QUOTIENT MINIMIZATION BY CONJUGATE GRADIENTS
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C    N   =  ORDER OF PROBLEM
C    X   =  INITIAL (APPROXIMATE?) EIGENVECTOR
C   EPS  =  MACHINE PRECISION
C&&& for Microsoft test replace with actual names
C   APR,BPR  ARE NAMES OF SUBROUTINES WHICH FORM THE PRODUCTS
C          V= A*X    VIA   CALL APR(N,X,V)
C          T= B*X    VIA   CALL BPR(N,X,T)
C   KPR  =  LIMIT ON THE NUMBER OF PRODUCTS (INPUT) (TAKES ROLE OF IPR)
C        =  PRODUCTS USED (OUTPUT)
C  Y,Z,T,G,A,B RE WORKING VECTORS IN AT LEAST N ELEMENTS
C   P0   =  APPROXIMATE EIGENVALUE (OUTPUT)
C   IPR  =  PRINT CHANNEL   PRINTING IF IPR.GT.0
C  STEP 0
      INTEGER N,LP,IPR,ITN,I,LIM,COUNT
      REAL X(N),T(N),G(N),Y(N),Z(N),PN,A(N),B(N)
      REAL EPS,TOL,P0,PA,XAX,XBX,XAT,XBT,TAT,TBT,W,K,D,V,GG,BETA,TABT,U
C  IBM VALUE   -  APPROX. LARGEST NUMBER REPRESENTABLE.
C&&&       PA=R1MACH(2)
      PA=1E+35
      LIM=KPR
      KPR=0
      TOL=N*N*EPS*EPS
C  STEP 1
  10  KPR=KPR+1
      IF(KPR.GT.LIM)RETURN
C  FIND LIMIT IN ORIGINAL PROGRAMS
      CALL APR(N,X,A)
      CALL BPR(N,X,B)
C  STEP 2
      XAX=0.0
      XBX=0.0
      DO 25 I=1,N
        XAX=XAX+X(I)*A(I)
        XBX=XBX+X(I)*B(I)
  25  CONTINUE
C  STEP 3
      IF(XBX.LT.TOL)STOP
C  STEP 4
      P0=XAX/XBX
      IF(P0.GE.PA)RETURN
      IF(IPR.GT.0)WRITE(IPR,963)KPR,P0
 963  FORMAT( 1H ,I4,' PRODUCTS, EST. EIGENVALUE=',1PE16.8)
C  STEP 5
      PA=P0
C  STEP 6
      GG=0.0
      DO 65 I=1,N
        G(I)=2.0*(A(I)-P0*B(I))/XBX
        GG=GG+G(I)**2
  65  CONTINUE
C  STEP 7
      IF(IPR.GT.0)WRITE(IPR,964)GG
 964  FORMAT(' GRADIENT NORM SQUARED=',1PE16.8)
      IF(GG.LT.TOL)RETURN
C  STEP 8
      DO 85 I=1,N
        T(I)=-G(I)
  85  CONTINUE
C  STEP 9
      DO 240 ITN=1,N
C  STEP 10
        KPR=KPR+1
        IF(KPR.GT.LIM)RETURN
        CALL APR(N,T,Y)
        CALL BPR(N,T,Z)
C  STEP 11
        TAT=0.0
        TBT=0.0
        XAT=0.0
        XBT=0.0
        DO 115 I=1,N
        TAT=TAT+T(I)*Y(I)
         XAT=XAT+X(I)*Y(I)
          TBT=TBT+T(I)*Z(I)
         XBT=XBT+X(I)*Z(I)
 115    CONTINUE
C  STEP 12
        U=TAT*XBT-XAT*TBT
        V=TAT*XBX-XAX*TBT
        W=XAT*XBX-XAX*XBT
        D=V*V-4.0*U*W
C  STEP 13
        IF(D.LT.0)STOP
C  MAY NOT WISH TO STOP
C  STEP 14
        D=SQRT(D)
        IF(V.GT.0.0)GOTO 145
        K=0.5*(D-V)/U
        GOTO 150
 145    K=-2.0*W/(D+V)
 150    COUNT=0
C  STEP 15
        XAX=0.0
        XBX=0.0
        DO 155 I=1,N
          A(I)=A(I)+K*Y(I)
          B(I)=B(I)+K*Z(I)
          W=X(I)
          X(I)=W+K*T(I)
          IF(W.EQ.X(I))COUNT=COUNT+1
          XAX=XAX+X(I)*A(I)
          XBX=XBX+X(I)*B(I)
 155    CONTINUE
C  STEP 16
        IF(XBX.LT.TOL)STOP
        PN=XAX/XBX
C  STEP 17
        IF(COUNT.LT.N)GOTO 180
        IF(ITN.EQ.1)RETURN
        GOTO 10
C  STEP 18
 180    IF(PN.LT.P0)GOTO 190
        IF(ITN.EQ.1)RETURN
        GOTO 10
C  STEP 19
 190    P0=PN
        GG=0.0
        DO 195 I=1,N
        G(I)=2.0*(A(I)-PN*B(I))/XBX
        GG=GG+G(I)**2
 195    CONTINUE
C  STEP 20
        IF(GG.LT.TOL)GOTO 10
C  STEP 21
        XBT=0.0
        DO 215 I=1,N
         XBT=XBT+X(I)*Z(I)
 215    CONTINUE
C  STEP 22
        TABT=0.0
        BETA=0.0
        DO 225 I=1,N
          W=Y(I)-PN*Z(I)
          TABT=TABT+T(I)*W
          BETA=BETA+G(I)*(W-G(I)*XBT)
 225    CONTINUE
C  STEP 23
        BETA=BETA/TABT
        DO 235 I=1,N
          T(I)=BETA*T(I)-G(I)
 235    CONTINUE
C  STEP 24
 240  CONTINUE
C  STEP 25
      GOTO 10
C  NO STEP 26  -  HAVE USED RETURN INSTEAD
      END
