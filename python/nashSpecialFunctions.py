#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 00:12:27 2021

@author: pcolsen
"""


def banana2(Bvec, a=1, b=100):
    """Calculate Two-dimensional Rosenbrock's Banana Function.

    This function implements Rosenbrock's Banana Function
    in two dimensions.  This function is a standard test
    function for minimization algorithms. It has a global
    minimum of 0 at the point (1, 1).  The classical starting
    point is (-1.2,1)


    Args:
        Bvec: a numeric sequence of length 2
        s, b: scale parameters, defauls 1, 100



    Returns:
        Real

    Ref:
        "An Automatic Method for finding the Greatest
        of Least Value of a Function" H. H. Rosenbrock
        The Computer Journal, Volume 3, Issue 3, 1960, Pages 175–184
    """

    if len(Bvec) != 2:
        raise TypeError("banana2: argument must be length 2")

    x = Bvec[0]
    y = Bvec[1]
    result = (a-x)**2 + b*(y-x**2)**2
    return result

def banana2n(Bvec, a=1, b=100):
    """Calculate even-dimensional Rosenbrock's Banana Function.

    This function implements Rosenbrock's Banana Function
    in even dimensions, 2*n.  This function is a standard test
    function for minimization algorithms. It has a global
    minimum of 0 at the point (1, 1, ..., 1).

    Args:
        Bvec: a numeric sequence of even length
        Classical starting point is a iteration
        of -1.2,1

    Returns:
        a number

    This is a sum of n two-dimensional banana functions.
    If banana2 works, this should as well

    """

    if (len(Bvec)%2) != 0:
        raise TypeError("banana2n: Bvec must have even length.")

    result = sum([banana2([Bvec[i], Bvec[i+1]], a, b)
                  for i in range(length(Bvec)/2)])
    return result


def bananaxn(Bvec):
    """Calculate more complex version of Rosenbrock's banana.

    This is a more involved variant of the banana function.
    This variant has exactly one minimum for n = 3, (at (1, 1, 1))
    and exactly two minima for 4<n<7, with the global minimum
    near (-1,1,1,...,1)

    Args:
        Bvec: a numeric vector

    Returns:
        a number
    """

    result = 100 * sum([((Bvec[i+1] - Bvec[i]**2)**2)
                        for i in range(len(Bvec) - 1)]
                       + [(1 + Bvec[i])**2
                          for i in range(len(Bvec) - 1)])
    return result


i

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 03:19:10 2021

@author: pcolsen
"""

from math import prod
from math import cos
from math import sqrt
import matplotlib.pyplot as plt

def griewank(Bvec):
    """Calculate Greiwank's function.

    Griewank's function is used to test
    maximization and minimization routines.

    Args:
        A numeric sequence

    Returns:
        Real
    """

    # if ((Bvec is not list)
    #     and (Bvec is not tuple)
    #     and (Bvec is not range)):
    #         raise TypeError("griewank: Bvec must be numeric sequence.")

    result = (1.0 + sum([x**2 for x in Bvec])/4000.0
              - prod([cos(Bvec[i]/sqrt(i+1)) for i in range(len(Bvec))]))
    return result



    return

if __name__ == "__main__":

    def testGriewank():
        """Test griewank in two dimensions.
        """

        plt.plot(range(-500, 500), [griewank([x]) for x in range(-500, 500)])

    testGriewank()






