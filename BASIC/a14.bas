3000 PRINT "SMEV.JC ALG 14 DEC 7 78"
3005 REM DIM A(N,N),V(N,N)
3010 FOR I=1 TO N
3015   FOR J=1 TO N
3020     LET V(I,J)=0
3025   NEXT J
3030   LET V(I,I)=1
3035 NEXT I
3040 LET I9=0
3045 LET I9=I9+1
3050 IF I9>30 THEN GOTO 3240
3055 LET N8=0
3060 FOR I=1 TO N-1
3065   FOR J=I+1 TO N
3070     LET P=A(I,J)+A(J,I)
3075     LET Q=A(I,I)-A(J,J)
3080     LET T=SQRT(P*P+Q*Q)
3085     IF T=0 THEN GOTO 3145
3090     IF Q<0 THEN GOTO 3120
3095     IF ABS(A(I,I))<ABS(A(I,I))+50*ABS(P) THEN GOTO 3105
3100     IF ABS(A(J,J))=ABS(A(J,J))+50*ABS(P) THEN GOTO 3145
3105     LET C=SQRT((T+Q)/(2*T))
3110     LET S=.5*P/(T*C)
3115     GOTO 3140
3120     LET S=SQRT((T-Q)/(2*T))
3125     IF P>0 THEN GOTO 3135
3130     LET S=-S
3135     LET C=.5*P/(T*S)
3140     IF 1<1+ABS(S) THEN GOTO 3155
3145     LET N8=N8+1
3150     GOTO 3220
3155     FOR K=1 TO N
3160       LET Q=A(I,K)
3165       LET A(I,K)=C*Q+S*A(J,K)
3170       LET A(J,K)=-S*Q+C*A(J,K)
3175     NEXT K
3180     FOR K=1 TO N
3185       LET Q=A(K,I)
3190       LET A(K,I)=C*Q+S*A(K,J)
3195       LET A(K,J)=-S*Q+C*A(K,J)
3200       LET Q=V(K,I)
3205       LET V(K,I)=C*Q+S*V(K,J)
3210       LET V(K,J)=-S*Q+C*V(K,J)
3215     NEXT K
3220   NEXT J
3225 NEXT I
3230 PRINT N8," SMALL ROTNS -- SWEEP",I9
3235 IF N8<N*(N-1)/2 THEN GOTO 3045
3240 PRINT "CONVERGED"
3245 FOR I=1 TO N
3250   LET D(I)=A(I,I)
3255 NEXT I
3260 RETURN 
