#!/bin/python3
import numpy

def FrankMat(n):
    Amat = numpy.empty(shape=(n,n), dtype='object')
    for i in range(1,n+1):
        print("i=",i)
        for j in range(1,i+1):
            print(j)
            Amat[i-1,j-1]=j
            Amat[j-1,i-1]=j
    return(Amat)



AA = FrankMat(4)
print(AA)
