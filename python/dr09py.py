# -*- coding: utf-8 -*-
"""
CNM Algorithm 09 test

J C Nash 2021-1-12
"""

import numpy
import math
import sys
def brspdmi(Avec, n):
# =============================================================================
# Bauer Reinsch inverse of symmetric positive definite matrix stored
#    as a vector that has the lower triangle of the matrix in row order
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

def FrankMat(n):
    Amat = numpy.array([ [ 0 ] * n ] * n) # numpy.empty(shape=(n,n), dtype='object')
    for i in range(1,n+1):
#        print("i=",i)
        for j in range(1,i+1):
#            print(j)
            Amat[i-1,j-1]=j
            Amat[j-1,i-1]=j
    return(Amat)

def smat2vec(Amat):
    n=len(Amat[0])
    n2=int(n*(n+1)/2)
    svec = [ None ] * n2
    k = 0
    for i in range(1,n+1):
        for j in range(1,i+1):
            svec[k]=Amat[i-1, j-1]
            k=k+1            
    return(svec)

def svec2mat(svec):
    n2=len(svec)
    n=int((-1+math.sqrt(1+8*n2))/2)
    print("matrix is of size ",n)
    Amat = numpy.array([ [ None ] * n ] * n)
    k = 0
    for i in range(1,n+1):
        for j in range(1,i+1):
            Amat[i-1, j-1] = svec[k]
            Amat[j-1, i-1] = svec[k]
            k=k+1            
    return(Amat)

# Main program
AA = FrankMat(4)
print(AA)
avec = smat2vec(AA)
print(avec)
n=len(AA[0])
vinv = brspdmi(avec, n)
## Computed inverse
##    2.00000 
##   -1.00000    2.00000 
##    0.00000   -1.00000    2.00000 
##    0.00000    0.00000   -1.00000    1.00000 



print(vinv)
Ainv = svec2mat(vinv)
print(Ainv)
print(AA)
print(numpy.dot(Ainv, AA))
