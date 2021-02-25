from calceps import *
import nashConstants

def hjmin(fminfn, B, intol=calceps()):

    '''Calculate the minimum of an n-dimentional function.
    =====================================================
    This function finds the minimum of an n-dimensional 
    target function using Hooke-Jeeves pattern search. To 
    do this it starts with some initial point in n-space,
    call it X, then looks for the minimum of the function on
    the vertices on the star defined by "X[i] + step" 
    centered on X.

    If there is a minimum on this start it then calculates 
    the vector, call it V that points from the original 
    point through the minimum vertex, takes a step along
    that vector in repeats the procedure centered at the new
    point. 

    If the procedure doesn't find a minimum on the star
    defined with the +step vertices, it tries the same 
    exploration on the -step verticies.  If this doesn't
    find a minimum, the it's in the vicinity of the
    minimum
    
    fminfn : the n-dimensional target function
    B : the n-dimensional point from which to begin the search
    intol : the tolerance used to find "close enough" to min.
    Defaults to calceps()

    returns a tuple (fmin, (tuple of argsmin))
    '''

    ########## Inner functions #############################
    # These functions are just to hide some complexity
    # Inner function to pretty print a vector
    def printVec(V):
        """Print elements of vec 7 per line
        """
        for v in V:
            index += 1
            print("{0} ".format(v))
            if 0 == (1index % 7):
                print("\n")
        return
    # End printVec #####

    # This function is syntactic sugar for exception trapping for
    # cases in which the Pascal programs would return noncomp
    # as a return argument. Python is call-by-value so we
    # can't return noncomp and use the try: idiom instead
    def try_fminfn_or_big(B):
        '''Returns value of fminfn or nashConstants.big if exception.'''
        try:
            fval = fminfn(B)
        except:
            fval =  nashConstants.big
        finally:
            return fval
        pass
     ########## End Inner Functions ########################  

    # fminfn must be a function object so it can be called
    if not (fminfn is function):
        raise TypeError("hjmin: target function fminfn must be a function")
        # Funrion exits here

    # Setting initial parameters
    result = None # This should never be returned.
    if intol < 0.0:
        intol = calceps() # The smalled meaningful value
    ifn = 1
    # fail signals if the algorithm has failed to find a min.
    fail = False
    stepsize = 0.0
    for i in range(length(B)):
        if stepsize <  nashConstants.stepredn*abs(B[i]):
            stepsize = nashConstants.stepredn*abs(B[i])
    if stepsize == 0:
        stepsize = nashConstants.stepredn
    X = B.copy()

    # Now we test the target function, fminfn.
    # The function must be computable.  The original Pascal algorithm
    # had fminfn return a flag if the function could not be computed.
    # The flag can't be returned in Python because Python calls by
    # value. We could return None, but that complicates fminfn
    # because it has to catch its own errors.  We'll trap them all.
    # Here we execute fminf with initial data, then catch exceptions. That
    # should yield the same result and it makes fminfn simpler. 
    try: # Catch errors
        # matches the "if noncomp" statement in Pascal version.
        # This iw the only place we check for computability
        # and the type of the return value.
        fval = fminfn(B)
        if fval is not (int or real):
            fail = True 
            # Return must be numeric and not complex
            raise TypeError("hjmin: fminfh: target function fminfn "
                                + "must return real or int,\n"
                                + "but was type {0}.".format(type(fval)))
        except: TypeError
            raise
        except: 
            # All other exceptions in target function.
            fail = True
            print("Error in fminfn in hjmin: function not computable at initial point.") 
            print("Unexpected error in fminfn in hjmin:", sys.exc_info()[0])
            raise
        # The rest of the function is in this else clause.
        else: # Matches the else of "if noncomp" in Pascal version.
            # We get here only if the fval is computable, numeric, and not complex
            # Print function value and arguments
            print("Initial function value: {0}\n".format(fval))
            printVec(B)
            # set initial state variables.
            Fmin = fval
            fold = fval 
            
            # Here's where we find the corner the cube with the minimum value.
            # We index through each point, first checking for a decrease
            # at a step up on the star.  If that doesn't work we try a
            # step down.  If that doesn't work we stay unmoved.
            # There is undoubtedly a more Pythonic way to do this.
            while (stepsize > intol):
                for i in range(length(B)): # 
                    temp = B[i]
                    # Start with a step up
                    B[i] = temp + stepsize
                    # We know that fminfn is computable at the original B,
                    # so we try it at a step up
                    fval = try_fminfn_or_big(B) # eval fminfn in inner function
                    if fval < Fmin:
                        Fmin = fval
                    # This leaves B[i] one step up
                    # Finding a min at step-up means there won't be one at
                    # step-down unless we're near a local max
                    else: 
                        # No point in looking here if we found a min at step-up
                        B[i] = temp - stepsize
                        fval = try_fminfn_or_big(B) # 
                        if fval < Fmin:
                            Fmin = fval # This leaves B[1] one step down.
                        else:
                            B[i]=temp # There was no lower min.

                # If Fmin < fold, then we have a downhill direction,
                # so we step in that direction,
                if Fmin < fold:
                    # Seve current position in X, then step out B
                    # This is a pattern move.
                    temp = [2.0*B[i]-X[i] for i in range(length(B))]
                    X = B.copy()
                    B = temp.copy()
                fold = Fmin
            else:
                # Here Fmin !< fold so we're either stuck in the same
                # spot, in which case we reduce the stepsize and try again,
                # or we're not in the same spot, in which case we center
                # on the new point
                
                # This replaces the repeat/until structure in the
                # Pascal routine. It will be true only if all the
                # point elements are equal.
                samepoint = all([B[i]==X[i] for i in range(length(B))])
;
                # If we are at the same point, shrink the stepsize,
                # then print the current value and point, 
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
        result =(Fmin, tuple(B)))
    return result
            

        
    
                 
    
    
            

    ## Unfinished.
