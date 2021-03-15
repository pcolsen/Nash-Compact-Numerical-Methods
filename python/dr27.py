#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 00:12:27 2021

@author: pcolsen
"""

banana2(Bvec):

    """Two-dimensional Rosenbrock's Banana Function

    This function implements Rosenbrock's Banana Function
    in two dimensions.  This function is a standard test
    function for minimization algorithms. It has a global 
    minimum of 0 at the point (1, 1).

    Args:
        Bvec: a numeric vector of length 2
        Classical starting point (-1.2,1)

    Returns:
        Real

    Ref:
        "An Automatic Method for finding the Greatest
        of Least Value of a Function" H. H. Rosenbrock
        The Computer Journal, Volume 3, Issue 3, 1960, Pages 175–184
    """

    if length(Bvec) != 2:
        raise TypeError("banana2: argument must be length 2")

    x = Bvec[0]
    u = Bvec[1]
    result = (1-x)**2 + 100*(y-x**2)**2
    return result

banana2n(Bvec):

    """Even-dimensional Rosenbrock's Banana Function

    This function implements Rosenbrock's Banana Function
    in n dimensions where n is even.  This function 
    is a standard test function for minimization algorithms.
    There are two ways to construct an n-dimensional 
    banana function.  This is an extension to n-dimensions
    consructed by summing a sequence of 2-dmensional 
    functions.

    Args:
        Bvec: a numeric vector of even length
        Classical starting point (-1.2,1)

    Returns:
        Real

    """


    if (length(Bvec)%2) != 0:
        raise TypeError("banana2n: Bvec must have even length.")

    result = sum([(100*(Bvec[2*i-1]**2 - Bvec[2*i])**2)
                  + ((Bvec[2*i-1] - 1)**2) for i in length(Bvec)/2])
    return result

