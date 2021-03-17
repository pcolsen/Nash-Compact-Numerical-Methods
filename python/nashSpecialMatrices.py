#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 22:52:24 2021

This module provides special matrices.

These matrices can be divided into two groups

Some have known propoeties.  These are useful for testing
matrix algorithms.

Others merely offer a convenient way to obtain matrices
that might be of practical interest.

Ducmentation for the testing matrices and their uses
can be found in "Compact Numerical Methods for Digital
Computers" by John Nash, particularly Appendix 1.
For several of these matrices further documentation
is available on the web.

Documentation for the convenience matrices is contained
in the function docstrings and in the documentation for
nashLib.

Example:
    hessian_matrix(size) returns a Hessian matrix of
    size-by-size.

    random_matrix(size) returns a size-by-size matrix
    of elements uniformly distributed on the interval [0,1)

All of the functions in this module use the matrix_broadcast
from the nashMatrixTools module

@author: pcolsen
"""

from math import floor
from random import uniform, random, shuffle, choice
from  nashMatrixTools import matrix_broadcast

def hilbert_matrix(size):
    """Return a matrix of size

    This function returns a Hilbert matrix in which
    A[i,jj = 1/ (i+j), where indices i and j are zero-based.


    Parameters
    ----------
    size : integer
        gives the size of the returned matrix

    Returns
    -------
    Hilbert matrix as list of lists.

    """
    def hilbert_element(a, b):
        return 1/(a+b+1)
    return matrix_broadcast(hilbert_element, size, size)


def dingdong_matrix(size):
    """Return a dingdong mastrix of size.

    This function returns a dingdong matrix of size in which
    elements are defined as A[i,j]=0·5/(n-i-j+1·5).
    The dingdong matrix is defined exclusively in "Compact
    Numerical Methods for Digital Computers" by John Nash,
    Appendix 1.

    Parameters
    ----------
    size : integer
        size of the desired matrix.

    Returns
    -------
    dingdong matrix as a list of lists

    """
    def dingdong_element(i, j):
        return 0.5/(size-i-j+0.5)
    return matrix_broadcast(dingdong_element, size, size)


def moler_matrix(size):
    """Return a moler matrix.

    Returns a matrix of size with elements
    A[i,i] =i, A[i,j]=min(i,j) for zero-based i,j.

    Parameters
    ----------
    size : integer
        size of the matrix to be returned.

    Returns
    -------
    Moler matrix of size.

    """
    def moler_element(i, j):
        if i == j:
            result = i
        else:
            result = min(i, j)-1
        return result
    return matrix_broadcast(moler_element, size, size)


def frank_matrix(size):
    """Returns a Frank matrix.

    Returns a size-by-size Frank matrix with elements
    Ai[i,j] = min(i,j).  A Frank matrix is
    reasonably well-behaved.

    Parameters
    ----------
    size : int
        THe dimension of the returned Frank matrix.

    Returns
    -------
    A size-by-size Frank matrix.

    """
    def frank_element(i,j):
        return min(i,j)
    return matrix_broadcast(frank_element, size, size)

def bordered_matrix(size):
    def bordered_element(i, j):
        if i == j:
            result = 1
        elif j == size-1:
            result = 1.0/(2**i)
        elif i == size-1:
            result = 1.0/(2**j)
        else:
            result = 0
        return result
    return matrix_broadcast(bordered_element, size, size)


def diagonal_matrix(size):
    def diagonal_element(i, j):
        if i == j:
            result = i
        else:
            result = 0
        return result
    return matrix_broadcast(diagonal_element, size, size)


def wilkinson_plus_matrix(size):
    def wilk_plus_elt(i, j):
        if i == j:
            result = abs(floor(size/2.0)-i)
        elif (j == i+1) or (i == j+1):
            result = 1
        else:
            result = 0
        return result
    return matrix_broadcast(wilk_plus_elt, size, size)


def wilkinson_minus_matrix(size):
    def wilk_plus_elt(i, j):
        if i == j:
            result = -floor(size/2.0)-i
        elif (j == i+1) or (i == j+1):
            result = 1
        else:
            result = 0
        return result
    return matrix_broadcast(wilk_plus_elt, size, size)


def ones_matrix(size):
    def ones_elt(i, j):
        return 1
    return matrix_broadcast(ones_elt, size, size)

def random_matrix(size):
    def uniform_elt(i, j):
        return random()
    return matrix_broadcast(uniform_elt, size, size)

def random_shuffled_matrix(size):
    base = list(range(size*size))
    shuffle(base)
    print(base)
    def shuffled_elt(i, j):
        return(base[i*size+j])
    return matrix_broadcast(shuffled_elt, size, size)


