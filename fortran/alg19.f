      SUBROUTINE A19NM(N,B,X,NX,NX2,NOCOM,IFN,VL,FUN,STEP,IPR)
C  ALGORITHM 19 NELDER-MEAD SIMPLEX FUNCTION MINIMIZATION
C  THIS VERSION MODIFIED 1989-04-25 IN ACCORD WITH ERRORS DISCOVERED SINCE
C  THE 1979 VERSION WAS RELEASED
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  N= NO. OF PARAMETERS
C  B= VECTOR CONTAINING STARTING POINT & RETURNING MINIMUM
C  X= WORKING ARRAY NX BY NX2
C  NX & NX2 ARE DIMENSIONS OF X.  NX.GE.(N+1), NX2.GE.(N+2)
C  NOCOM=LOGICAL FLAG SET .TRUE. IF INITIAL POINT INFEASIBLE OR STEP=0.0
C  IFN = LIMIT TO NO. OF FUNCTION EVALUATIONS (INPUT)
C      = FUNCTION EVALUATIONS USED (OUTPUT)
C  VL  = FUNCTION VALUE AT MINIMUM
C  FUN = NAME OF FUNCTION SUBROUTINE  P= FUN(N,B,NOCOM)
C  STEP= STEPSIZE FOR CONSTRUCTING INITIAL SIMPLEX
C  IPR = PRINTER CHANNEL.  PRINTING IF IPR.GT.0
C  STEP0
      LOGICAL NOCOM
      INTEGER N,NX,NX2,N1,N2,IFN,LIFN,H,L,C,IPR
C
      REAL B(N),X(NX,NX2),STEP,P,VL,SIZE,ALPHA,BETA,GAMMA,T,VH,SS,VNEXT
      C=N+2
      N1=N+1
      ALPHA=1.0
      BETA=0.5
      GAMMA=2.0
      LIFN=IFN
C  IBM VALUE
C&&&       BIG=R1MACH(2)
      BIG=1.0E+35
C  STEP 1
      NOCOM=.FALSE.
      P=FUN(N,B,NOCOM)
      IFN=1
C  PRECAUTION FOR NULL STEP
      IF(STEP.EQ.0.0)NOCOM=.TRUE.
      IF(NOCOM)RETURN
      IF(IPR.GT.0)WRITE(IPR,960)P
 960  FORMAT('0INITIAL FUNCTION VALUE=',1PE16.8)
C  STEP 2
      X(N1,1)=P
      DO 10 I=1,N
        X(I,1)=B(I)
  10  CONTINUE
      L=1
      SIZE=0.0
C  STEP 3
      DO 40 J=2,N1
C  STEP 4
        DO 20 I=1,N
          X(I,J)=B(I)
  20    CONTINUE
        T=STEP
C  STEP 5
  25    X(J-1,J)=B(J-1)+T
C  STEP 6
        IF(X(J-1,J).NE.B(J-1))GOTO 30
        T=10.0*T
        GOTO 25
C  STEP 7
  30    SIZE=SIZE+ABS(T)
C  STEP 8
  40  CONTINUE
C  STEP 9
  45  SS=SIZE
C  STEP 10
      DO 60 J=1,N1
C  STEP 11
        IF(J.EQ.L)GOTO 60
C  STEP 12
        DO 50 I=1,N
          B(I)=X(I,J)
  50    CONTINUE
      NOCOM=.FALSE.
        P=FUN(N,B,NOCOM)
        IF(NOCOM)P=BIG
        IFN=IFN+1
        X(N1,J)=P
C  STEP 13
  60  CONTINUE
C  STEP 14
  65  L=1
      H=1
C****       NEXT=1
      VL=X(N1,1)
      VH=VL
C  STEP 15
      DO 80 J=2,N1
C  STEP 16
        T=X(N1,J)
C  STEP 17
        IF(T.GE.VL)GOTO 70
        VL=T
        L=J
C  STEP 18
  70    IF(T.LT.VH)GOTO 80
C****         NEXT=H
        H=J
        VH=T
C  STEP 19
  80  CONTINUE
C  PRINTOUT
      IF(IPR.GT.0)WRITE(IPR,950)IFN,VL,VH
 950  FORMAT( 8H #EVALS=,I4,19H SIMPLEX LOW + HIGH,1P2E16.8)
      IF(IFN.GT.LIFN)GOTO 400
C  STEP 20
      IF(VL.EQ.VH)RETURN
C**** FOLLOWING STATEMENT IS NEW
      VNEXT=VL+BETA*(VH-VL)
C     THIS SETS THE FUNCTION VALUE AT THE 'NEXT TO HIGHEST POINT'
C     IN AN APPROXIMATE WAY WITHOUT THE NECESSITY OF SEARCHING
C     FOR THIS POINT ITSELF.
C  STEP 21
      DO 100 I=1,N
        T=-X(I,H)
        DO 90 J=1,N1
          T=T+X(I,J)
  90    CONTINUE
        X(I,C)=T/N
 100  CONTINUE
C  STEP 22
      DO 110 I=1,N
        B(I)=(1.0+ALPHA)*X(I,C)-ALPHA*X(I,H)
 110  CONTINUE
      NOCOM=.FALSE.
      P=FUN(N,B,NOCOM)
      IFN=IFN+1
      IF(NOCOM)P=BIG
C  STEP 23
      IF(P.LT.VL)GOTO 350
C  STEP 24
      IF(P.LT.VNEXT)GOTO 390
C  STEP 25
      IF(P.GE.VH)GOTO 270
C  STEP 26
      DO 265 I=1,N
        X(I,H)=B(I)
 265  CONTINUE
      X(N1,H)=P
C  STEP 27
 270  DO 275 I=1,N
        B(I)=(1.0-BETA)*X(I,H)+BETA*X(I,C)
 275  CONTINUE
      NOCOM=.FALSE.
      P=FUN(N,B,NOCOM)
      IFN=IFN+1
      IF(NOCOM)P=BIG
C  STEP 28
      IF(P.LT.X(N1,H))GOTO 390
C  STEP 29
      SIZE=0.0
      DO 320 J=1,N1
C  STEP 30
        IF(J.EQ.L)GOTO 320
C  STEP 31
        DO 310 I=1,N
          X(I,J)=BETA*(X(I,J)-X(I,L))+X(I,L)
          SIZE=SIZE+ABS(X(I,J)-X(I,L))
 310    CONTINUE
C  STEP 32
 320  CONTINUE
C  STEP 33
      IF(SIZE.LT.SS)GOTO 45
C  STEP 34
      GOTO 400
C  STEP 35
 350  DO 355 I=1,N
        T=GAMMA*B(I)+(1.0-GAMMA)*X(I,C)
        X(I,C)=B(I)
        B(I)=T
 355  CONTINUE
      X(N1,C)=P
C  STEP 36
      NOCOM=.FALSE.
      P=FUN(N,B,NOCOM)
      IFN=IFN+1
      IF(NOCOM)P=BIG
C  STEP 37
      IF(P.LT.X(N1,C))GOTO 390
C  STEP 38
      DO 385 I=1,N
        B(I)=X(I,C)
 385  CONTINUE
      P=X(N1,C)
C  STEP 39
 390  DO 395 I=1,N
        X(I,H)=B(I)
 395  CONTINUE
      X(N1,H)=P
      GOTO 65
C  STEP 40
 400  IF(L.EQ.0)NOCOM=.TRUE.
      IF(NOCOM)RETURN
      DO 410 I=1,N
        B(I)=X(I,L)
 410  CONTINUE
      RETURN
      END