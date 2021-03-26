/*     CHOLTEST.C
C  {Program to test the timings of the Choleski decompostion by
C   three variants of the same algorithm}
C   MAIN
*/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <conio.h>

typedef struct {
   double    *element;     /* the individual elements (element[index])*/
   int N;
   int tmsize; /* usually 1+N*(N+1)/2 */
} tmatrix;

/* the following routines create a matrix of the appropraiate size and     */
/* return a pointer to it.  */

tmatrix   *new_tmatrix(int N); /* create a new triangular mat  */

/* The following matrix routines do not create storage for the results     */
void        disposetmatrix(tmatrix *A);                 /* free the space    */

int nchol(tmatrix *A);
int pchol(tmatrix *A);
int hchol(tmatrix *A);
void dumptmat(tmatrix *A);

main ()
{
tmatrix *A;
tmatrix *B;
float TIMES[5][3];
clock_t TTOT, T0, T1;
int  IFL;
int IREP, IRR, N, NN, N2, I, J, K;
FILE *out;


  IREP = 100;      /* normally 100 ?? */
/*  C     10 REPETITIONS  */
  printf("\n CHOLTEST.FOR -- timing study of Choleski variants \n");
  for (NN=1; NN<=5; NN++)
  {   N=10*NN;
      printf("  TEST AT ORDER= %d \n",N); /* getch(); */
      N2=N*(N+1)/2;
/*      printf("about to allocate matrix B \n"); getch(); */
      B=new_tmatrix(N);
/*      printf("B allocated \n"); */
      if(!B) exit(999);
/* C    {Build Moler matrix} */
      for (I=1; I<=N; I++)  /* DO 10 */
      {
	for (J=1; J<=I; J++) /*  DO 8 J=1,I */
	{
	  K = I*(I-1)/2 + J;
/*	  printf("K=%d ",K,":");  getch();  */
	  B->element[K]=(double)(J-2);
/*	  printf("%g ",B->element[K]);  getch(); */
	}
/*	printf(" old K=%d ",K); getch(); */
	B->element[K]=(double)I;
/*	printf("%g ",B->element[K]);  getch(); */
      }
/*      getch(); printf("\n"); */
/*      printf("Temp halt \n"); getch();
      exit(22); */
/*      dumptmat(B);*/

/* C    {Matrix now built == begin timings} */
       A=new_tmatrix(N); /* create matrix */
/* C    {First Nash}  */
      printf(" NASH \n");
      TTOT=0.0;
      for (IRR=1; IRR<=IREP; IRR++) /*    DO 20 IRR=1,IREP */
      {
	for (K=1; K<=N2; K++) A->element[K] = B->element[K];
	T0 = clock();
	if (nchol(A) == 1)
	{printf("\n  ** ERROR ** failure in NCHOL.C  \n");
	 return 1;
	};
/*	printf("DECOMPOSED:\n"); */
/*	dumptmat(A);*/
	T1 = clock();
	TTOT=TTOT+(T1-T0);
      } /*  20   CONTINUE */
      TIMES[NN-1][0] = TTOT/(float)CLK_TCK;

      printf("  PRICE \n");
      TTOT=0;
      for (IRR=1; IRR<=IREP; IRR++) /* DO 30 IRR=1,IREP */
      {
	for (K=1; K<=N2; K++) A->element[K] = B->element[K];
	T0 = clock();
	if (pchol(A) == 1)
	{printf("\n  ** ERROR ** failure in PCHOL.C  \n");
	 return 2;
	};
	T1 = clock();
	TTOT=TTOT+(T1-T0);
      } /*  20   CONTINUE */
      TIMES[NN-1][1] = TTOT/(float)CLK_TCK;
      printf("  HEALY \n");

      TTOT=0;
      for (IRR=1; IRR<=IREP; IRR++) /* DO 30 IRR=1,IREP */
      {
	for (K=1; K<=N2; K++) A->element[K] = B->element[K];
	T0 = clock();
	if (hchol(A) == 1)
	{printf("\n  ** ERROR ** failure in HCHOL.C  \n");
	 return 2;
	};
	T1 = clock();
	TTOT=TTOT+(T1-T0);
      } /*  20   CONTINUE */
      TIMES[NN-1][2] = TTOT/(float)CLK_TCK;
      disposetmatrix(A);
      disposetmatrix(B);
/*  100   CONTINUE  */
   }
/* C END MAIN LOOP */
   printf("\n ORDER     NASH         PRICE      HEALY \n");
   for (NN=1; NN<=5; NN++)
   { N=10*NN;
/*      WRITE(OUT,905)N, TIMES(NN,1), TIMES(NN,2), TIMES(NN,3)
905   FORMAT(1H ,I6,4X,3F10.2)
120   CONTINUE                */
     printf("%d      %f     %f     %f \n",N,TIMES[NN-1][0],TIMES[NN-1][1],TIMES[NN-1][2]);
   }
   out = fopen("CHOLMC","wt");
   fprintf(out,"\n ORDER     NASH         PRICE      HEALY \n");
   for (NN=1; NN<=5; NN++)
   { N=10*NN;
     fprintf(out,"%d      %f     %f     %f \n",N,TIMES[NN-1][0],TIMES[NN-1][1],TIMES[NN-1][2]);
   }
  return 0;
} /* end main */


int nchol(tmatrix *A)
/* C NASH */
{
float S;
int I, J, J1, L, L2, LK, L2K, K, N;

{
  N = A->N;
/*
  printf("\n  begin NASH \n");
  for (I=1; I<=N; I++)
  { for (J=1; J<=I; J++)
    { printf("%f ",A->element[I*(I-1)/2+J]); }
      printf("\n");
  }
  getch(); */

for (J=1; J<=N; J++) /*   DO 20 J=1,N  */
  { L=J*(J+1)/2;
    if (J>1)    /*      IF(J.EQ.1) GO TO 11  */
    {  J1=J-1;
       for (I=J; I<=N; I++) /*   DO 10 I=J,N */
       { L2=I*(I-1)/2+J;
	 S=A->element[L2];
	 for (K=1; K<=J1; K++) /*  DO 5 K=1,J1 */
	 { L2K=L2-K;
	   LK=L-K;
	   S=S-A->element[L2K]*A->element[LK];
	 }   /* 5       CONTINUE */
	 A->element[L2]=S;
       } /* 10       CONTINUE */
/* 11       IF(A(L).LE.0.0) GO TO 21  /*
    }
    if (A(L)<=0.0)
    { return 1 }
    else
    {  S=sqrt(A->element[L]);
       for (I=J; I<=N; I++) /*    DO 19 I=J,N  */
       {
	L2=I*(I-1)/2+J;
	A->element[L2]=A->element[L2]/S;
       } /* 19       CONTINUE */
    } /* end if */
 } /* 20       CONTINUE */
/*  printf("\n  end NASH \n");
  for (I=1; I<=N; I++)
  { for (J=1; J<=I; J++)
    { printf("%f ",A->element[I*(I-1)/2+J]); }
      printf("\n");
  }
  getch(); */
 return 0;
}
} /* end nchol */

int pchol(tmatrix *A)
/* C PRICE */
{
float S, T;
int I, I1, I2, J, J1, J2, L, L2, LK, L2K, K, M, N;
  N = A->N;

/*  printf("\n  begin PRICE \n");
  for (I=1; I<=N; I++)
  { for (J=1; J<=I; J++)
    { printf("%f ",A->element[I*(I-1)/2+J]); }
      printf("\n");
  }
  */
   J2=0;
   for (J1=1; J1<=N; J1++)  /*    DO 20 J1=1,N  */
   {
       I2=J2;
       for (I1=J1; I1<=N; I1++) /*   DO 10 I1=J1,N  */
       {
	 L=I2+J1;
	 S=A->element[L];
	 J=1;
/* 4       IF(J.EQ.J1) GO TO 5  */
	 while (J<J1)
	 {
	  K=I2+J;
	  M=J2+J;
	  S=S-A->element[K]*A->element[M];
	  J=J+1;
	 }
/*       GOTO 4  */
/* 5       IF(I1.EQ.J1) GO TO 6 */
	 if (I1 > J1)
	 {
	   A->element[L]=S/T;
	 }
/*       GO TO 7  */
	 else
	 {
/* 6       IF(S.LE.0.0) GO TO 21 */
	   if (S <= 0.0)  return 1;
	   T=sqrt(S);
	   A->element[L]=T;
	  } /* end else */
/* 7       I2=I2+I1 */
	 I2=I2+I1;
/* 10       CONTINUE  */
       } /* end inner for */
       J2=J2+J1;
/* 20       CONTINUE */
    } /* end outer for */
/*  printf("\n  end PRICE \n");
  for (I=1; I<=N; I++)
  { for (J=1; J<=I; J++)
    { printf("%f ",A->element[I*(I-1)/2+J]); }
      printf("\n");
  }
  getch(); */
    return 0;
} /* end pchol */


int hchol(tmatrix *A)
/* C HEALY */
{
float S,W;
int I, ICOL, IROW, J, J1, L, L2, LK, L2K, K, M, N;
  N = A->N;
/*  printf("\n  begin HEALY \n");
  for (I=1; I<=N; I++)
  { for (J=1; J<=I; J++)
    { printf("%f ",A->element[I*(I-1)/2+J]); }
      printf("\n");
  }  */

  J=1;
  K=0;
/*  DO 10 ICOL=1,N */
  for (ICOL=1; ICOL<=N; ICOL++)
  {
    L=0;
/*  DO 11 IROW=1,ICOL  */
    for (IROW=1; IROW<=ICOL; IROW++)
    {
      K=K+1;
      W=A->element[K];
      M=J;
/*      DO 12 I=1,IROW  */
      for (I=1; I<=IROW; I++)
      {
       L=L+1;
/*       IF(I.EQ.IROW) GO TO 13  */
       if (I<IROW)
       {
	W=W-A->element[L]*A->element[M];
	M=M+1;
       }
       else break;
      } /* end for I */
/* 12   CONTINUE  */
/* 13   IF(IROW.EQ.ICOL) GO TO 14 */
      if (IROW != ICOL)
      {
      if(A->element[L]<=0.0) return 1;
/*         IF(A(L).LE.0.0) GO TO 21  */
      A->element[K]=W/A->element[L];
/* 11   CONTINUE  */
      } /*  end if */
    }  /* end for IROW */
/* 14   IF(W.LE.0.0) GO TO 21  */
    if(W <= 0.0) return 1;
    A->element[K]=sqrt(W);
    J=J+ICOL;
/*  10   CONTINUE  */
  } /* end for ICOL */
/*  printf("\n  end HEALY \n");
  for (I=1; I<=N; I++)
  { for (J=1; J<=I; J++)
    { printf("%f ",A->element[I*(I-1)/2+J]); }
      printf("\n");
  }
  getch(); */
   return 0;
}

/* allocates the storage for the matrix of the specified size */
tmatrix   *new_tmatrix(int N)
{
int i;
tmatrix   *mytmat;
/*   printf("new_tmatrix: allocate matrix of order %d \n",N); getch();*/
   if (!N) return NULL; /* only matrices with dimensions */
   mytmat=(tmatrix *)malloc(sizeof(tmatrix));
   if (!mytmat) return NULL;
/*   printf("new_tmatrix: malloc OK \n"); getch();*/
   mytmat->tmsize = 1 + N*(N+1) / 2; /* not using row zero */
/*   printf("new_tmatrix: tmsize = %d \n",mytmat->tmsize); getch(); */
     /* fill out with appropriate # of elements */
      mytmat->element=(double *)calloc(mytmat->tmsize,sizeof(double));
      if (!mytmat->element) {
	 free((void *)mytmat);
         return NULL;
      }
   mytmat->N = N;
/*   printf(" mytmat->N = %d \n",mytmat->N); getch(); */
   return mytmat;
}  /* end new_tmatrix */

/* free all the space taken by the matrix */
void  disposetmatrix(tmatrix *A)
{
int i;
   if (!A) return;
      free((void *)A->element);   /* free the elements           */
      free((void *)A);            /* the structure overhead  */
} /* end disposetmatrix */

/* void dumptmat(tmatrix *A)
{
int i,j;
    if (!A){
	printf("Dumptmat: matrix NOT defined \n");
	exit(1);
	}
    for (i=1; i<=A->N; i++)
    {
      for (j=1; j<=i; j++) printf("%g ",A->element[i*(i-1)/2+j]);
      printf("\n");
    }
    return;
}*/ /* end dumptmat */
