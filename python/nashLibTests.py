from nashMatrixTools import *
from nashSpecialMatrices import *
from nashLib import *

def check_nashGivens():
    print("Checking the nashGivens algorithm.")
    
    testMatrix = [[1, 1, 1, 1, 1],
                [1, 2, 2, 2, 2],
                [1, 2, 3, 3, 3],
                [1, 2, 3, 4, 4],
                [1, 2, 3, 4, 5]]
    print("Test Matrix\n")
    pprint_matrix(testMatrix)
    print("\n\n")
    
    Q,R = givens(testMatrix)
    
    print("Calculated Q matrix (rounded)\n")
    pprint_matrix(round_matrix(Q))
    print("\n")
    
    check_Q = [[0.447, -0.894, 0.000,  0.000,  0.000],
                   [0.447,   0.224,   -0.866,  0.000,   0.000],
                   [0.447,   0.224,  0.289,   -0.816,  0.000],
                   [0.447,   0.224,   0.289,   0.408,   -0.707],
                   [0.447,   0.224,   0.289,   0.408,   0.707]]
        
    print("Known-correct Q matrix (roiunded)\n")
    pprint_matrix(check_Q)
    print("\n\n")
    
    print("Calculated R matrix (rounded)\n")
    pprint_matrix(R)
    print("\n")
    
    check_R = [[2.236,   4.025,   5.367,   6.261,   6.708],
                   [0.000,   0.894,   1.565,   2.012,   2.236],
                   [0.000,   0.000,   0.866,   1.443,   1.732],
                   [0.000,   0.000,   0.000,   0.816,   1.225],
                   [0.000,   0.000,   0.000,   0.000,   0.707]]
        
    print("Known-correct R matrix (rounded)\n")
    pprint_matrix(R)
    print("\n\n")
    
    print("The next matrix should bew an identity.\n")
    pprint_matrix(round_matrix(multiply_matrices(Q, transpose_matrix(Q))))
    
    return


