C&&& A23
C  TEST ALG. 23 J C NASH AUG 1978
C  J.C. NASH   JULY 1978, APRIL 1989
      LOGICAL NOCOM
      INTEGER M,N,N2,LIM,IG,NIN,NOUT
      REAL PM,TOL,A(10),C(10),B(4),X(4),V(4),D(4),F(7)
      EXTERNAL WRES,WDRES
      N=4
      M=7
C  I/O CHANNELS
      NIN=5
      NOUT=6
  10  READ(NIN,900)LIM,TOL
 900  FORMAT(I5,F15.8)
      WRITE(NOUT,950)LIM,TOL
 950  FORMAT(' PROBLEM WOOD4  FN EVAL LIMIT=',I5,'  LAMBDA TOL=',E16.8)
      IF(LIM.LE.0)STOP
      N2=N*(N+1)/2
      READ(NIN,901)B
 901  FORMAT(4F10.5)
      WRITE(NOUT,953)B
 953  FORMAT(1H ,4F15.5)
      CALL A23MRT(N,B,M,TOL,A,C,N2,X,V,D,WRES,WDRES,NOCOM,PM,LIM,IG,
     #F,0)
      IF(NOCOM)WRITE(NOUT,951)
 951  FORMAT(' FAILURE')
      WRITE(NOUT,952)PM,LIM,IG
 952  FORMAT(' CONV. TO',1PE16.8,' IN',I5,' FNS &',I4,' DERS')
      WRITE(NOUT,953)B
      WRITE(NOUT, 902)
 902  FORMAT(' ')
      GOTO 10
      END
      SUBROUTINE A23MRT(N,B,M,TOL,A,C,N2,X,V,D,RES,DRES,NOCOM,P0,IFN,
     #IG,F,IPR)
C  ALGORITHM 23  MODIFIED MARQUARDT NONLINEAR SUM OF SQUARES
C    MINIMISATION
C  J.C. NASH   JULY 1978, FEBRUARY 1980, APRIL 1989
C     N  =  NO. OF PARAMETERS TO BE ADJUSTED
C     B  =  INITIAL POINT (SET OF PARAMETERS
C     M  =  NO. OF RESIDUALS
C    TOL =  RESET VALUE FOR MARQUARDT PARAMETER LAMBDA
C    A.C =  WORKING VECTORS OF N2 ELEMENTS
C  X,V,D =  WORKING VECTORS OF N  ELEMENTS
C  RES   =  NAME OF FUNCTION TO CALCULATE RESIDUAL NO. I
C           RVAL=RES(N,B,I,NOCOM)
C  DRES  =  NAME OF SUBROUTINE TO CALCULATE DERIVATIVES OF RESIDUAL I
C           CALL DRES(N,B,I,D)
C  NOCOM =  LOGICAL FLAG SET .TRUE. IF INITIAL POINT INFEASIBLE
C  P0    =  MINIMAL VALUE OF SUM OF SQUARES (OUTPUT)
C  IFN  =  LIMIT ON FUNCTION EVALUATIONS (INPUT)  (SUM OF SQUARES)
C        =  COUNT OF FUNCTION EVALUATIONS (OUTPUT)
C  IG    =  COUNT OF DERIVATIVE EVALUATIONS
C  F     =  WORKING VECTOR OF LENGTH M USED TO SAVE RESIDUALS
C  IPR  =  PRINT CHANNEL  IPR.GT.0 FOR PRINTING.
C  STEP 0
      LOGICAL NOCOM
      INTEGER N,M,N2,IFN,IG,LIM,I,J,Q,IJ,J1,COUNT
      REAL B(N),X(N),V(N),D(N),A(N2),C(N2),F(M)
      REAL  S,TOL,INC,DEC,LAMBDA,PHI,P,P0
C  FOR SAFETY RESET N2
      N2=N*(N+1)/2
C  PHI - NASH ADDITION TO MARQUARDT ALGORITHM
      PHI=1.0
C  INCREASE AND DECREASE FACTORS
      INC=10.0
      DEC=0.4
      LIM=IFN
      IFN=0
      IG=0
      LAMBDA=TOL
C  STEP 1
      P=0.0
C  BETTER DONE DOUBLE PRECISION
      IFN=IFN+1
      DO 15 I=1,M
        F(I)=RES(N,B,I,NOCOM)
        IF(NOCOM)RETURN
        P=P+F(I)**2
  15  CONTINUE
C  STEP 2
  20  IG=IG+1
      LAMBDA=LAMBDA*DEC
      P0=P
      IF(IPR.GT.0)WRITE(IPR,959)IG,IFN,P0
 959  FORMAT( 6H ITN #,I4, 8H EVALN *,I4,13H  SUMSQUARES=,1PE16.8)
C  STEP 3
      DO 34 J=1,N2
        A(J)=0.0
  34  CONTINUE
      DO 36 J=1,N
        V(J)=0.0
  36  CONTINUE
C  STEP 4
      DO 48 I=1,M
        CALL DRES(N,B,I,D)
        S=F(I)
        DO 46 J=1,N
           V(J)=V(J)+S*D(J)
           Q=J*(J-1)/2
           DO 44 K=1,J
            IJ=Q+K
             A(IJ)=A(IJ)+D(J)*D(K)
  44       CONTINUE
  46    CONTINUE
  48  CONTINUE
C  STEP 5
      DO 54 J=1,N2
        C(J)=A(J)
  54  CONTINUE
      DO 56 J=1,N
        D(J)=B(J)
  56  CONTINUE
C  STEP 6
  60  DO 68 J=1,N
        Q=J*(J+1)/2
        A(Q)=C(Q)*(1.0+LAMBDA)+PHI*LAMBDA
        X(J)=-V(J)
        IF(J.EQ.1)GOTO 68
        J1=J-1
        DO 66 I=1,J1
          IJ=Q-I
          A(IJ)=C(IJ)
  66    CONTINUE
  68  CONTINUE
C  STEP 7
      NOCOM=.FALSE.
      CALL A7CH(A,N2,N,NOCOM)
      IF(NOCOM)GOTO 130
C  STEP 8
      CALL A8CS(A,N2,X,N)
C  STEP 9
      COUNT=0
      DO 95 I=1,N
        B(I)=D(I)+X(I)
        IF(B(I).EQ.D(I))COUNT=COUNT+1
  95  CONTINUE
C  STEP 10
      IF(COUNT.EQ.N)RETURN
C  STEP 11
      IFN=IFN+1
      IF(IFN.GT.LIM)GOTO 140
      NOCOM=.FALSE.
      P=0.0
      DO 115 I=1,M
        F(I)=RES(N,B,I,NOCOM)
        IF(NOCOM)GOTO 130
        P=P+F(I)**2
 115  CONTINUE
C  STEP 12
      IF(P.LT.P0)GOTO 20
C  STEP 13
 130  LAMBDA=LAMBDA*INC
      IF(LAMBDA.EQ.0.0)LAMBDA=TOL
      GOTO 60
C  RESET PARAMETERS
 140  DO 144 I=1,N
        B(I)=D(I)
 144  CONTINUE
      RETURN
      END
      FUNCTION WRES(N,B,I,NOCOM)
C  J.C. NASH   JULY 1978, APRIL 1989
      LOGICAL NOCOM
      INTEGER N,I
      REAL B(N),FV,D(4)
      NOCOM=.FALSE.
      CALL WOOD4(FV,D,B,I,.FALSE.)
      WRES=FV
      RETURN
      END
      SUBROUTINE WDRES(N,B,I,D)
C  J.C. NASH   JULY 1978, APRIL 1989
      INTEGER N,I
      REAL B(N),D(N),FV
      CALL WOOD4(FV,D,B,I,.TRUE.)
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
      SUBROUTINE WOOD4(FVAL,D,X,I,MODE)
C  J.C. NASH   JULY 1978, APRIL 1989
C  WOOD'S 4 PARAMETER FUNCTION
C  FVAL  =  FUNCTION VALUE OR RESIDUAL VALUE AT POINT I
C   D    =  DERIVATIVES OF RESIDUAL I
C   X    =  POINT
C   I    =  OBSERVATION NO.  IF 0  THEN  COMPUTE SUM OF SQUARES OR GRAD
C  MODE  =  F  FN OR  RES, T  DERIVS. OR GRADIENT
      LOGICAL MODE
      INTEGER I
      REAL X(4),D(4)
C  HILLSTROM EXPRESSIONS
      IF(MODE) GOTO 500
      IF(I.GT.0)GOTO 250
C   SUM OF SQUARES TOTAL FN
      FVAL=100.0*(X(2)-X(1)**2)**2 + (1.0-X(1))**2
      FVAL=FVAL + 90.0*(X(4)-X(3)**2)**2
      FVAL=FVAL + (1.0-X(3))**2
      FVAL=FVAL + 10.1*((X(2)-1.0)**2 + (X(4)-1.0)**2)
      FVAL=FVAL + 19.8*(X(2)-1.0)*(X(4)-1.0)
      RETURN
C   RESIDUALS
 250  GOTO (310,320,330,340,350,360,370),I
 310  FVAL=10.0*(X(2)-X(1)**2)
      RETURN
 320  FVAL=1.0-X(1)
      RETURN
 330  FVAL=SQRT(90.0)*(X(4)-X(3)**2)
      RETURN
 340  FVAL=1.0-X(3)
      RETURN
 350  FVAL=SQRT(0.2)*(X(2)-1.0)
      RETURN
 360  FVAL=SQRT(0.2)*(X(4)-1.0)
      RETURN
 370  FVAL=SQRT(9.9)*(X(2)+X(4)-2.0)
      RETURN
C    DERIVATIVES
 500  IF(I.GT.0)GOTO 750
C   GRADIENT OF FN
      D(1)=-400.0*(X(2)-X(1)**2)*X(1) - 2.0*(1.0-X(1))
      D(2)= 200.0*(X(2)-X(1)**2) + 20.2*(X(2)-1.0) + 19.8*(X(4)-1.0)
      D(3)=-360.0*X(3)*(X(4)-X(3)**2) - 2.0*(1.0-X(3))
      D(4)= 180.0*(X(4)-X(3)**2) + 20.2*(X(4)-1.0) + 19.8*(X(2)-1.0)
      RETURN
C   RESIDUAL DERIVS AT OBSN I
 750  GOTO (810,820,830,840,850,860,870), I
 810  D(1)= -20.0*X(1)
      D(2)= 10.0
      D(3)=  0.0
      D(4)=  0.0
      RETURN
 820  D(1)= -1.0
      D(2)=  0.0
      D(3)=  0.0
      D(4)=  0.0
      RETURN
 830  D(1)=  0.0
      D(2)=  0.0
      D(3)= -2.0*SQRT(90.0)*X(3)
      D(4)=  SQRT(90.0)
      RETURN
 840  D(1)= 0.0
      D(2)=  0.0
      D(3)= -1.0
      D(4)=  0.0
      RETURN
 850  D(1)=  0.0
      D(2)= SQRT(0.2)
      D(3)=  0.0
      D(4)=  0.0
      RETURN
 860  D(1)=  0.0
      D(2)=  0.0
      D(3)=  0.0
      D(4)= SQRT(0.2)
      RETURN
 870  D(1)=  0.0
      D(2)= SQRT(9.9)
      D(3)=  0.0
      D(4)= SQRT(9.9)
      RETURN
      END
C&&&   300  1.0
C&&&   -3.0      -1.0      -3.0      -1.0
C&&&   300  0.1
C&&&   -3.0      -1.0      -3.0      -1.0
C&&&   300  0.0001
C&&&   -3.0      -1.0      -3.0      -1.0
C&&&   300  0.0000001
C&&&   -3.0      -1.0      -3.0      -1.0
C&&&    10 0.0001
C&&&   -3.0      -1.0      -3.0      -1.0
C&&&     0
