##### Unfinished

def hjmin(targetf, B, intol):

    '''Calculate the minimum of an n-dimentional function.
    =====================================================
    This function finds the minimum of an n-dimensional 
    target function using Hooke-Jeeves pattern search. To 
    do this it starts with some initial point in n-space,
    call it X, then explores the vertices of an n-cube
    centered on X. It then calculates the vector, call it V,
    that points from the vertex with the maximum value of
    the target function, Xmax, to the one with the minimum
    value, Xmin.  It then conducts a one-dimensional search
    along that extended vector, finds the minimum, and starts
    again.
    
    targetf : the n-dimensional target function
    B : the n-dimensional point from which to begin the search
    intol : the tolerance used to find "close enough" to min.
    '''

    if not (targetf is function):
        raise TypeError("hjmin: targetf must be a function")
    
    if intol < 0.0:
        intol = calceps()
    # fail signals if the algorithm has failed to find a min.
    fail = False
    # the step used to define the n-cube
    stepsize = min([b*nashConstants.stepredn*abs(b) for b in B])
    # but stepsize must not be zero
    stepsize = max(0,nashConstants.stepredn)
    
            

    ## Unfinished.
