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






