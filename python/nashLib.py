from nashSpecialMatrices import *
from nashMatrixTools import *
from math import sqrt


def givens(A):
    # Python indexes from 0, so rows are indexed 0 to nRow-1 
    nRow = len(A)
    # Python indexes from 0, so rows are indexed 0 to nRow-1 
    nCol = len(A[0])
    mn = min(nRow,nCol) # Step 0
    # mn is the minimum of nRow and nCo1 and gives the maximum size
    # of the triangular matrix resulting from the reduction. Note 
    # that the decomposition is still valid when nRow<nCol, but 
    # the R matrix is then trapezoidal, i.e. an upper trianglular 
    # matrix with added columns on the right.
    eps = calceps()
    Q = identity_matrix(mn) # We'll calculate the rotations in Q
    # Step 1: main loop on diagonals of triangle
    for j in range(mn): # looping over rows with j
        # Step 2 loop on column elements with i
        for k in range(j+1, nRow):
            c = A[j][j]
            s = A[k][j]
            b = max(abs(c), abs(s))
            # Step 4 normalise elements to avoid over- or under-flow
            if b > 0: # 
                c = c/b
                s = s/b
                p = sqrt(c*c+s*s)
                s = s/p
                if abs(s) >= eps: # Step 5
                    c = c/p;
                    # Step 7 -- rotation of A
                    for i in range(nCol):
                        p = A[j][i]
                        A[j][i] = c*p+s*A[k][i]
                        A[k][i] = -s*p+c*A[k][i]
                    # Step 8- rotation of Q. Note: nRow not nCo1.
                    for i in range(nRow):
                        p = Q[i][j]
                        Q[i][j] = c*p+s*Q[i][k]
                        Q[i][k] = -s*p+c*Q[i][k]
    return [Q, A]

print("working_givens.py has been loaded.")
 
