#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 18:16:19 2021.

@author: pcolsen
"""

from itertools import tee

def dot(x, y):
    """Return the dot product of two sequences."""
    if len(x) != len(y):
        raise TypeError("dot: vectors were different lengths.")
    result = sum(map(lambda x, y: x*y, x, y))
    return result


def vector_scale(vector, scaler):
    """Scale a vector."""
    return list(map(lambda x: scaler*x, vector))

def scale_matrix(matrix, scaler):
    """Scale a matrix or a vector."""
    if matrixp(matrix):
        result = [vector_scale(v, scaler) for v in matrix]
    else:
        result = vector_scale(matrix, scaler)
    return result



def flip_matrix(Vector):
    """Zips a sequence or sequence of sequences.

    flip takes a list-of-lists [1,1,1,...,1],
    [2,2,2,....,2],...,[n,n,n,...,n] and zips them to
    be [1,2,3,...n], 1,2,3,...,n]....

    This function is used to prepare the right-hand
    matrix for matrix multiplication.  It doesn't
    transpose a matrix, it just re-slices it vertically
    instead of horizontally

    Args:
    -----
        Vector: a swquence.

    Returns
    -------
        A list of lists zipped
    """
    return list(map(lambda x: list(x), list(zip(*Vector))))


def matmult(M1, M2):
    """Multiply and MxN matrix by NxP.

    There's a special case when M2 is a
    vector, so P=1.

    This function multiplies two matrices by
    using the 'flip' function to re-ravel the
    right-hand matrix, then taking the outer-product
    of the rows and columns

    Args:
    -----
        M1: MxN matrix
        M2: NxP matrix

    Returns
    -------
        MxP matrix product
    """

    if len(M1[0]) != len(M2):
        raise TypeError("malmult: inconformable matrices.")
    if matrixp(M2):
        result = outerproduct(dot, M1, flip(M2))
    else:
        result = list(map(lambda v,s : vector_scale(v,s), M1, M2))
    return result

    """Apply func to the outer prod of rows in M1 to cols in M2.

    Note that M1 and M2 do not generally need to be conformable
    for matrix multiplication.  That is required only in the
    case when the rows of M1 and the columns of M2 are
    arguments to a function that needs equal-length
    sequences as arguments, such as dot.

    args:
        func: a function of two vectors

    return [[func(u,v) for v in flip(M2)] for u in M1]
    """
    return [[func(u,v) for v in M2] for u in M1]


def vectorp(thing):
    """Return True if thing is a list or a tuple."""
    return type(thing) in (list, tuple)


def matrixp(thing):
    """Return True if Matrix, False if Vector, None if neither."""
    # Must be a seqnence to be either
    if not vectorp(thing):
        result = None
    else: # A natrix needs a sequence first element
        result = vectorp(thing[0])
    return result


def ravel_vector(vector, rowlength):
    """Turn a vector into a matrix with rowlength.


    Parameters
    ----------
    vector : a vector
        the source vector for the matrix.
    rowlength : scaler
        the length of the rows in the resultant matrix

    Returns
    -------
        a matrix with rows of rowlength
    """
    if (len(vector)%rowlength) != 0:
        raise ValueError("ravel: rowlength not a divisor of vector length.")
    result = [[vector[i] for i in list(range(j, j+rowlength))]
              for j in list(range(0, len(vector), rowlength))]
    return result


def unravel_matrix(matrix):
    """Unravel a matrix into a vector.


    Parameters
    ----------
    matrix : a list-of-lists matrix
        the mastrix to be unraveled into a vector

    Returns
    -------
    A vector of the contents of the unraveled matrix

    """
    return [inner for outer in matrix for inner in outer]

def matrix_broadcast(func, rows, columns):
    """Construct a matrix by applying func to indices.

    This routine constructs a matrix by applying the
    function func to the indices i,j at every position
    in the matrix.

    func can be any function that returns a value in
    response to two arguments.  The value need not depend
    explicitly on i or j---for exapmple it may be a constant,
    and it need not even be a number---for example it could
    be a string, or even a vector, which would yield a
    tensor of order three.

    Examples of how this function is used can be found
    in the file nashSpecialMatrices.py

    Parameters
    ----------
    func : a function of two arguments
        the result of this function is the value applied to
        the element at the argument indices, e.g. "i+j"
    rows: int
        The number of rows in the result matrix.
        Must be > 0
    columns: int
        The number of columns in the result matrix.
        Must be > 0.

    Returns
    -------
    A list-of-lists matrix with the elements determined by func

    """
    return [[func(i,j) for j in range(columns)] for i in range(rows)]


def negate_matrix(mat):
    return scale_matrix(mat, -1)

def join_matrices(M1, M2):
    """Join two matrices row-by-row.

    Takes two matrices and adds the elements of each row of
    the right-hand matrix onto the end of each row of the left.
    The matrices must have the same number of rows.  They may
    have different numbers of columns

    Example:
        m1 = [[1, 2], [5, 6]]
        m2 = [[3, 4], [7, 8]]

        join_matrices(m1, m2)
        [[1, 2, 3, 4],
         [5, 6, 7, 8]]

    Parameters
    ----------
    M1 : list-of-lists matrix
    M2 : list-of-lists matrix

    Raises
    ------
    ValueError
        The two matrices have the same number of rows.

    Returns
    -------
    M11 : list-of-lists matrix
        Returns the two matrices joined row-by-row.
    """
    if len(M1) != len(M2):
        raise ValueError("joing_matrices: args must have same number of rows.")
    else:
        newmat = copy_matrix(M1)
        [newmat[i].extend(M2[i]) for i in range(len(newmat))]
    return newmat


def stack_matrices(M1, M2):
    """Join two matrices column-by-column.

    Takes two matrices and adds the elements of each column of
    the right-hand matrix onto the bottom of each column of the left.
    The matrices must have the same number of columns.  They may
    have different numbers of rows.

    Example:
        m1 = [[1,3,5,7]]
        m2 = [[2,4,6,8]]

    stack_matrices(m1,m2)
        [[1, 3, 5, 7],
         [2, 4, 6, 8]]

    Parameters
    ----------
    M1 : list-of-lists matrix
    M2 : list-of-lists matrix

    Raises
    ------
    ValueError
        The two matrices have the same number of columns.

    Returns
    -------
    M11 : list-of-lists matrix
        Returns the two matrices joined row-by-row.
    """
    if len(M1[0]) != len(M2[0]):
        raise ValueError("join_matrices: args must have same number of columns.")
    else:
        M11 = M1.copy()
        M11.extend(M2)
    return M11



def quad_tile_matrices(M1, M2, M3, M4):
    """Return list-of-lists 2x2 matrix tiling 4 args.

    This function takes the four square argument matrices and
    produces another square lists-of-lists matrix
    composed by tiling them.

    Example
    -------
    quad_tile_matrices(M1, M2, M3, M4) yields

    [M1 M2
     M3 M4]

    m1 = [[1,2],[5,6]]
    m2 = [[3,4],[7,8]]
    m3 = [[9,10],[13,14]]
    m4 = [[11,12],[15,16]]

    quad_tile_matrices(m1,m2,m3,m4) yields

    [[1, 2, 3, 4],
     [5, 6, 7, 8],
     [9, 10, 11, 12],
     [13, 14, 15, 16]]


    Parameters
    ----------
    M1, M2, M3, M4: list-of-lists square matrices

    Returns
    -------
    list-of-lists matrix
        Returns the list-of-lists matrix created by
        tiling the four square argument matrices.

    """

    return stack_matrices(join_matrices(M1, M2),
                          join_matrices(M3, M4))

def copy_matrix(M1):
    """Return a new copy of the argument matrix.

    The function returns a copy of a list-of-lists
    matrix.

    There are two reasons for this function.

    First: Direct assignment, A = B, returns
    a new pointer to the old matrix, not a new
    matrix.  Changeing A[i][j] also changes
    B[i][j].  The test A is B will return true.

    Second, using the copy methos, A.copy()
    renerats a new instance of the ouuter list
    of matrix A, but does not generate new copies
    of its content.  A = B will return False, but
    A[i] B[i] will return True.  That means that
    while operations on the outer lists will
    be change only one matrix, operations on the
    original inner lists (rows) will effect both
    vectors. This can cause horribly confusing
    errors.


    Parameters
    ----------
    M1 : list-of-lists matrix
        This is the matrix to be copied.

    Returns
    -------
    list-of-lists matrix
        This is a separate copy of the original matrix.
        M1 is copy_matrix(M1) will return False, and
        M1[i] is copy(M1)[i] will return False as well.

    """
    return [M1[i][:] for i in range(len(M1))]


def transpose_matrix(M):
    """Transpose list-of-lists matrix.

    Parameters
    ----------
    M : list-of-lists matrix
        M is the matrix to be transposed.

    Returns
    -------
    list-of-lists matrix
        Returns the transpose of the argument list-of-lists matrix.

    """
    return [[M[j][i] for j in range(len(M))] for i in range(len(M[1]))]









