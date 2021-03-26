/*
 *	svd()
 *
 *	Singular Value Decomposition using orthogonalization by plane rotations.
 *	Algorithm #1 from the book by J. C. Nash, "Compact Numerical Methods
 *	for Computers: Linear Algebra and Function Minimisation," Wiley 1979.
 *
 *	This version performs the orthogonalization by rows.
 *
 *	This routine determines the singular value decomposition A = U S V'
 *	of a real m by n rectangular matrix.
 *
 * Usage:	int svd(), rank, m, n;
 *		double **u, **s, **v;
 *
 *		rank = svd( u, s, v, m, n);
 *
 *     On input:
 *
 *        u must be either an m-by-m identity matrix, or NULL if computation
 *	  of u is not desired.
 *
 *	  s is a 1-by-n matrix of zeros.
 *
 *        v contains the m-by-n input matrix to be decomposed.
 *
 *        m is the number of rows of v and s, and the order of u.
 *
 *        n is the number of columns of v and s.
 *
 *     On output:
 *
 *	  The return value is rank(a).
 *
 *        s contains the singular values in decreasing order.
 *
 *        If u was not NULL on input, u contains the matrix u (m-by-m
 *	  orthogonal) of the decomposition.
 *
 *        v contains the matrix v (n-by-n orthogonal) of the decomposition.
 *
 *     Note: if m > n, then on output, only the first n rows of v
 *	     are significant.
 */

#include <stdio.h>	/* for NULL */
#include <math.h>
#include "mat.h"

#define Q(i,j) q[i][j]
#define U(i,j) u[i][j]
#define V(i,j) v[i][j]
#define Z(i)   z[0][i]

static double ip(double *x, double *y, int n) /* inner product of two vectors */
{
    double r = 0.0;

    while( n--)
	r += *x++ * *y++;

    return r;
}

int svd( double *u[], double *z[], double *v[], int m, int n)
{
    int i, j, k, count, vrows;
    double rtol, ztol, p, q, r, c, s, d;
    double *tmp;

    vrows = m < n ? m : n;

    rtol = ztol = 10 * n * geteps();	/* tolerance for rank */
    ztol *= ztol;			/* tolerance for zero */

    if( m > n) {			/* first perform QR if m > n */
	qr( u, v, m, n);
	if( u)
	    mattrans1( u, m);
    }

    for( i = 0; i < vrows; ++i)
	Z(i) = ip( v[i], v[i], n);	/* store squared row norms in z */

    do {
	count = (vrows*(vrows-1)) / 2;

	for( j = 0; j < vrows-1; ++j)
	  for( k = j+1; k < vrows; ++k) { /* make rows j & k of v orthogonal */

	    q = Z(j); r = Z(k);

	    if( q < r) {				/* reverse order */

		tmp = v[j]; v[j] = v[k]; v[k] = tmp;
		if( u) {
		    tmp = u[j]; u[j] = u[k]; u[k] = tmp;
		}
		Z(j) = r; Z(k) = q;

	    } else {
		    p = ip( v[j], v[k], n);

		    if(  q*r < ztol || (p*p)/(q*r) < ztol)	/* orthogonal */
			--count;
		    else {					/* rotate */
			    q -= r;
			    d = sqrt(4.0*p*p+q*q);
			    c = sqrt((d+q)/(2.0*d));
			    s = p/(d*c);
			    for( i = 0; i < n; ++i) {
				p = V(j,i);
				V(j,i) = c * p + s * V(k,i);
				V(k,i) = c * V(k,i) - s * p;
			    }
			    if( u) for( i = 0; i < m; ++i) {
				p = U(j,i);
				U(j,i) = c * p + s * U(k,i);
				U(k,i) = c * U(k,i) - s * p;
			    }
			    Z(j) = ip( v[j], v[j], n);
			    Z(k) = ip( v[k], v[k], n);
		    }
	    }
	}
    } while( count);

    for( count = 0, i = 0; i < vrows; ++i) {	/* normalize rows of v */
	Z(i) = sqrt(Z(i));
	if( Z(i) >= ztol && Z(i)/Z(0) >= rtol) {
	    ++count;
	    for( j = 0; j < n; ++j)
		V(i,j) /= Z(i);
	}
    }

    if( m > 1 && u)			/* transpose u and v */
	mattrans1( u, m);

    mattrans1( v, n);

    if( count < n) {

	if( count == 0) {	/* v should be identity */

	    zeros( v, n, n);
	    eye( v, n, n);

	} else {		/* use QR to get remaining columns of v */

	    double **q, **r;

	    if( (q = matcreate( n, n)) == NULL)
		error("matcreate() failed in svd()","");

	    if( (r = matcreate( n, count)) == NULL)
		error("matcreate() failed in svd()","");

	    eye( q, n, n);
	    matcopy( r, v, n, count);
	    qr( q, r, n, count);

	    for( i = 0; i < n; ++i)
		for( j = count; j < n; ++j)
		    V(i,j) = Q(i,j);

	    matfree( q); matfree( r);
	}
    }

    return count;
}
