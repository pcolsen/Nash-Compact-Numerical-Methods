##### Unfinished

def hjmin(fminfn, B, intol):

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
    # Inner function to pretty print a vector
    def printVec(V):
        """Print elements of vec 7 per line
        """
        for v in V:
            index += 1
            print("{0} ".format(v))
            if 0 == (1index % 7):
                print("\n")
    # End printVec
       

    # fminfn must be a function object so it can be called
    if not (fminfn is function):
        raise TypeError("hjmin: targetf must be a function")
    
    if intol < 0.0:
        intol = calceps() # The smalled meaningful value
    ifn = 1
    # fail signals if the algorithm has failed to find a min.
    fail = False
    # the step used to define the n-cube
    stepsize = min([b*nashConstants.stepredn*abs(b) for b in B if b > 0])
    # but stepsize must not be zero
    stepsize = max(0,nashConstants.stepredn)
    # X will hold current B values while calculating search directions.
    X = B.copy()
    # The function must be computable.  The original Pascal algorithm
    # had fminfn return a flag if the function could not be computed.
    # The flag can't be returned in Python because Python calls by
    # value. We could return None, but that complicates this function.
    # Here we execute fminf with initial data, then catch exceptions. That
    # should yield the same result and it makes fminfn simpler. 
    try: # Catch errors if 
        fval = fminfn(B)
    except ValueError as err:
        print("Error in fminfn in hjmin: function not computable at initial point.") 
        print("ValueError in fminfn in hjmin: {0}".format(err))
        raise
    except TypeError as err:
        print("Error in fminfn in hjmin: function not computable at initial point.") 
        print("TypeError in fminfn in hjmin: {0}".format(err))
        raise
    except:
        print("Error in fminfn in hjmin: function not computable at initial point.") 
        print("Unexpected error in fminfn in hjmin:", sys.exc_info()[0])
        raise

    # Function could be computed
    print("Initial function value: {0}".format(fval))
    index = 0
    # Print elements of B when called.
    printVec(B)
    Fmin = fval # fval is fminfn at B when called.
    fold = fval # Keeps track the minimum fminfn at each cycle.
    # Here's where we find the corner the cube with the minimum value.
    # There is undoubtedly a more Pythonic way to do this.
    while (stepsize > intol):
           for i in range(length(B)): # 
               temp = B[i]
               # Start with a step up
               B[i] = temp + stepsize
               # We know that fminfn is computable at the original B,
               # so we try it at a step up
               try:
                   fval = fminfn(B)
               except: # We won't handle this explicitly, just plug a max value 
                   fval = nashConstants.big
               if fval < Fmin:
                   Fmin = fval
                   # This leaves B[i] one step up
                   # Finding a min at step-up means there won't be one at
                   # step-down unless we're near a local max
                else: 
                    # No point in looking here if we found a min at step-up
                    B[i] = temp - stepsize
                    try:
                        ifn += 1
                        fval = fminfn(B)
                    except:
                        fval = nashConstants.big
                    if fval < Fmin:
                        Fmin = fval # This leaves B[1] one step down.
                    else:
                        B[i]=temp # We never got a lower min.
                pass
            pass # end of for loop
            if Fmin < fold: # We've found a downhill direction.
                for i in range(length(b)):
                    # B is the new downhill vector, X is the last position
                    temp = 2.0*B[i]-X[i] # Why the factor of 2?
                    X[i] := B[i]
                    B[i] := temp
                fold = Fmin
            else
                samepoint = True
                i = 1
                while (samepoint and i<=n):
                    if B[i] != X[i]:
                        samepoint = False;
                        i += 1;
                if samepoint:
                    stepsize = stepsize*nashConstants.stepredn
                    print("hjmin: stepsize now: {0}".format(stepsize)
                              + ",  Best fn value: {0}".format(fmin)
                              + ",  after: {0}".format(ifn)
                              + " function evals.\n\n")
                    # Print B using inner function.
                    printVec(B)
                else:
                    B = X.copy()
    
    
        print("Converged to Fmin= {0}".format(Fmin)
                  ", after {0}".format(ifn)+" evaluations");
            
  return (Fmin, X)          
            

        
    
                 
    
    
            

    ## Unfinished.
