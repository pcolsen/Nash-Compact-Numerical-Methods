# -*- coding: utf-8 -*-
def brspdmi(Avec, n):
# =============================================================================
# Bauer Reinsch inverse of symmetric positive definite matrix stored
#    as a vector that has the lower triangle of the matrix in row order
#  alg09.py
# =============================================================================
    print(Avec)
    X = numpy.array([ 0 ] * n) # zero vector x
    for k in range(n, 0, -1):
        s  =  Avec[0];
        #print("s=",s)
        if (s > 0.0) :
            m  =  1;
            for i in range(2,n+1):
                q = m
                m = m+i
                t = Avec[q]
                X[i-1] = -t/s
                if i>k :
                    X[i-1] = -X[i-1]
           #     print("i, q, m:", i, q, m)
                for j in range((q+2), m+1):
           #         print(j)
           #         print("j-q-1=",j-q-1)
           #         print(X[j-q-1])
                    Avec[j-i-1] = Avec[j-1]+t*X[j-q-1]
                q = q-1
                Avec[m-1] = 1.0/s
            for i in range(2, n+1):
                print("i ",i)
                Avec[q+i-1] = X[i-1]
        else :
            print("Matrix is singular")
            sys.exit()
        print(k,":",Avec)
    return(Avec)
