C&&& A18
C  TEST ALG 18  A18RF USING FORSYTHE,MALCOLM + MOLER
C  J.C. NASH   JULY 1978, APRIL 1989
      EXTERNAL FUNS
      INTEGER IPR,COUNT,NBIS,NIN
      LOGICAL NOCOM
      REAL U,V,TOL
      NIN=5
      IPR=6
   1  READ(NIN,901)COUNT,NBIS,U,V,TOL
 901  FORMAT(2I5,3F10.5)
      WRITE(IPR,910)COUNT,NBIS,U,V,TOL
 910  FORMAT(' TEST- COUNT=',I5,' NBIS=',I5,' U=',F15.5,' V=',F15.5,
     *' TOL=',F15.10)
      IF(COUNT.LE.0)STOP
      CALL A18RF(U,V,COUNT,NBIS,FUNS,TOL,NOCOM,IPR)
      IF(NOCOM)WRITE(IPR,961)
 961  FORMAT(' FAILURE')
      WRITE(IPR,962)U,V,COUNT
 962  FORMAT(' ROOT U=',1PE16.8,'  F(U)=',E16.8,' AFTER',I4,' EVALNS')
      GOTO 1
      END
      REAL FUNCTION FUNS(X,NOCOM)
C  J.C. NASH   JULY 1978, APRIL 1989
      LOGICAL NOCOM
      REAL X
      NOCOM=.FALSE.
      IF(X.LT.-5.0.OR.X.GT.10.0)NOCOM=.TRUE.
C  DEFINE RANGE OF INTEREST
      IF(NOCOM)RETURN
      FUNS=X*(X*X-2.0)-5.0
      RETURN
      END
      SUBROUTINE A18RF(U,V,COUNT,NBIS,FUN,TOL,NOCOM,IPR)
C  ALGORITHM 18
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C  BISECTION  FALSE-POSITION ROOTFINDER
C  U,V DEFINE INTERVAL OF INTEREST
C  ON OUTPUT ROOT IS AT U,  FUN(U)=V
C  COUNT = LIMIT ON FUNCTION EVALUATIONS (INPUT)
C        = NO. OF EVALUATIONS USED (OUTPUT)
C  NBIS = INTERVAL BETWEEN BISECTIONS - BISECTION PERFORMED
C         EVERY NBIS STEPS
C  FUN = NAME OF FUNCTION
C     CALLING SEQUENCE    FVAL=FUN(X,NOCOM)
C  NOCOM=LOGICAL FLAG SET .TRUE. IF ALGORITHM CANNOT PROCEED
C  TOL = CONVERGENCE TOLERANCE ON LENGTH OF INTERVAL - STEP 9A
C  IPR = PRINTER CHANNEL -- NO OUTPUT UNLESS IPR.GT.0
      LOGICAL NOCOM
      INTEGER COUNT,LIMIT,NBIS,IPR,NBC
      REAL U,V,TOL,FU,FV,B,FB
C  STEP 0
      LIMIT=COUNT
      COUNT=0
      NOCOM=.FALSE.
      COUNT=COUNT+1
      FU=FUN(U,NOCOM)
      IF(NOCOM)RETURN
      FV=FUN(V,NOCOM)
      COUNT=COUNT+1
      IF(NOCOM)RETURN
C  STEP 1
      IF(FU*FV.GT.0.0)NOCOM=.TRUE.
C  ALTERNATIVE USES SIGN FUNCTION TO AVOID UNDERFLOW
      IF(NOCOM)RETURN
C  STEP 2
  20  B=(U*FV-V*FU)/(FV-FU)
C  STEP 3
  30  IF(B.GT.U)GOTO 40
      B=U
      V=FU
      RETURN
  40  IF(B.LT.V)GOTO 50
C  STEP 4
C  QUESTION .LE..GE. STEPS 3 + 4 OF MS. SEEMS WRONG
      U=V
      V=FV
      RETURN
C  STEP 5
  50  IF(IPR.GT.0)WRITE(IPR,970)COUNT,U,FU,V,FV
      COUNT=COUNT+1
 970  FORMAT(1H ,I4,11H EVALNS, F(,1PE16.8,2H)=,E16.8,4H  F(,E16.8,
     *2H)=,E16.8)
      IF(COUNT.GT.LIMIT)GOTO 110
      FB=FUN(B,NOCOM)
      IF(NOCOM)RETURN
C  STEP 6
      IF(FB*FU.GT.0.0)GOTO 80
C  STEP 7
      FV=FB
      V=B
      GOTO 90
C  STEP 8
  80  FU=FB
      U=B
C  STEP 9A
  90  IF((V-U).GT.TOL)GOTO 95
      V=FU
      RETURN
 95   NBC=(COUNT-2)/NBIS
      IF(NBIS*NBC.NE.(COUNT-2))GOTO 20
C  STEP 10
      IF(IPR.GT.0)WRITE(IPR,971)COUNT
 971  FORMAT(' BISECTION AT EVALN #',I5)
      B=U+(V-U)*0.5
      GOTO 30
C  NOTE ACTION IN CASE OF FUNCTION EVALUATION OVER-RUN
 110  NOCOM=.TRUE.
      COUNT=LIMIT
      V=FU
      RETURN
      END
C&&&     5    5   0.0       3.0     0.000001
C&&&    40    5   0.0       3.0     0.0
C&&&    80    1   0.0       3.0     0.0
C&&&    40    5   0.0       3.0     0.001
C&&&    40    5   0.0       1.0     0.001
C&&&     0
