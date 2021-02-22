#######  Unfinished

import random

def constantVector(num, const=1):
    '''Returns a num-long vector of const
    ======================================
    num = length of the desired vector (>0)
    const = constant value of all elements.
    

    Most users could do this perfectly well by
    themselves.  This is syntactic sugar to
    make the operation explicit.'''

    if not num is integer:
        raise TypeError("Array length, num, must be an integer")
    if not num > 0:
        raise ValueError("Array length, num, must be > 0")

    vec = num * [const]
    return vec

def randomVector(num, upper=1, lower=0):
    '''Returns a n-long random vector in [0,bound)
    ==============================================
    num = length of the desired vector (>0)
    upper = upper bound on random numbers
    lower = lower bound on random numbers
    returned elements will be in [lower, upper)

    Most users can do this for themselves.  
    This is syntactic sugar to make the operation
    explicit'''

    randvec = [random.uniform(lower, upper) for x in range(num)]
    return randvec

def vectorin(num = None);
    '''Lets user enter a vector
    ============================
    num = length of input vector

    This function is designed to replace
    the vectorin procedure in the original 
    Nash-lib code.  There are much nicer to 
    do these things, but

    
