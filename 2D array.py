# -*- coding: utf-8 -*-
from numpy import asarray
# TODO: Replace all TODO comments (yes, this one too!)
ENERGY_LEVEL = [100, 113, 110,
                85, 105, 102, 86,
                63, 81, 101, 94,
                106, 101, 79,
                94, 90, 71]
# ====================================================
# The brute force method to solve first problem


def find_significant_energy_increase_brute(A):
    """
    Return a tuple (i,j)
    where A[i:j] is the most
    significant energy increase
period.
    time complexity = O(n^2)
    """
    # TODO
    max = 0
    initial = 0
    final = 0
    for x in range(len(A)):
        for y in range(x + 1, len(A)):
            if max < A[y] - A[x]:
                max = A[y] - A[x]
                initial = x
                final = y
    return (initial, final)


print('Brute method')
z = find_significant_energy_increase_brute(ENERGY_LEVEL)
print(z)

# ==========================================
# The recursive method to solve first problem


def find_significant_energy_increase_recursive(A):
    """
    Return a tuple (i,j) where A[i:j] is
     the most significant energy increase
period.
    time complexity = O (n logn)
    """
    # TODO
    b = max_difference(A, 0, len(A) - 1)
    return (b[1], b[2])


def max_difference(x, initial, final):
    if (initial == final) | (initial == final - 1):
        return[x[final] - x[initial], initial, final]

    else:
        middle = (initial + final) // 2
        a = max_difference(x, initial, middle)
        b = max_difference(x, middle + 1, final)
        c = max_difference_crossover(x, initial, middle, final)
        if a[0] > b[0]:
            if a[0] > c[0]:
                return a
            else:
                return c
        else:
            if b[0] > c[0]:
                return b
            else:
                return c


def max_difference_crossover(x, initial, middle, final):
    d = 0
    lmax = -999999
    rmax = -999999
    lpos = initial
    rpos = final
    for i in range(initial, middle):
        d = x[middle] - x[i]
        if d > lmax:
            lmax = d
            lpos = i
    for i in range(middle + 1, final):
        d = x[i] - x[middle + 1]
        if d > rmax:
            rmax = d
            rpos = i
    return [lmax + rmax, lpos, rpos]


print('Recursive method')
z = find_significant_energy_increase_recursive(ENERGY_LEVEL)
print(z)

# =========================================
# The iterative method to solve first problem


def find_significant_energy_increase_iterative(A):
    """
    Return a tuple (i,j) where A[i:j] is
    the most significant energy increase
period.
    time complexity = O(n)
    """
    # TODO

    low = A[0]
    max = A[1] - A[0] - 1
    initial = 0
    final = 0
    for a in range(1, len(A)):
        if A[a] - low > max:
            max = A[a] - low
            initial = A.index(low)
            final = a
        if A[a] < low:
            low = A[a]
    return (initial, final)


print('Iterative method')
z = find_significant_energy_increase_iterative(ENERGY_LEVEL)
print(z)

# =============================================
# The Strassen Algorithm to do the matrix multiplication


def sumup(x, y):
    n = len(x)
    result = [[0 for i in range(0, n)] for j in range(0, n)]
    for i in range(0, n):
        for j in range(0, n):
            result[i][j] = x[i][j] + y[i][j]
    return result

    # subtracts two matrices


def difference(x, y):
    n = len(x)
    result = [[0 for i in range(0, n)] for j in range(0, n)]
    for i in range(0, n):
        for j in range(0, n):
            result[i][j] = x[i][j] - y[i][j]
    return result


def square_matrix_multiply_strassens(x, y):
    """
    Return the product AB of matrix multiplication.
    Assume len(A) is a power of 2
    """

    x = asarray(x)

    y = asarray(y)

    assert x.shape == y.shape

    assert x.shape == x.T.shape

    assert (len(x) & (len(x)-1)) == 0, "A is not a power of 2"
    n = len(x)

    if n == 1:
        C = [[0 for j in range(0, n)] for i in range(0, n)]
        for i in range(0, n):
            for j in range(0, n):
                C[i][j] = x[i][j] * y[i][j]
        return C
    else:  # dividing the input matrices A and B
        new_n = int(n / 2)

    a11 = [[0 for i in range(0, new_n)] for j
           in range(0, new_n)]
    a12 = [[0 for i in range(0, new_n)] for j
           in range(0, new_n)]
    a21 = [[0 for i in range(0, new_n)]
           for j in range(0, new_n)]
    a22 = [[0 for i in range(0, new_n)]
           for j in range(0, new_n)]

    b11 = [[0 for i in range(0, new_n)]
           for j in range(0, new_n)]
    b12 = [[0 for i in range(0, new_n)]
           for j in range(0, new_n)]
    b21 = [[0 for i in range(0, new_n)]
           for j in range(0, new_n)]
    b22 = [[0 for i in range(0, new_n)]
           for j in range(0, new_n)]

    aTemp = [[0 for i in range(0, new_n)]
             for j in range(0, new_n)]
    bTemp = [[0 for i in range(0, new_n)]
             for j in range(0, new_n)]

    for i in range(0, new_n):
        for j in range(0, new_n):
            a11[i][j] = x[i][j]
            a12[i][j] = x[i][j + new_n]
            a21[i][j] = x[i + new_n][j]
            a22[i][j] = x[i + new_n][j +
                                     new_n]

            b11[i][j] = y[i][j]
            b12[i][j] = y[i][j + new_n]
            b21[i][j] = y[i + new_n][j]
            b22[i][j] = y[i + new_n][j +
                                     new_n]

    aTemp = sumup(a11, a22)
    bTemp = sumup(b11, b22)
    p1 = square_matrix_multiply_strassens(aTemp, bTemp)

    aTemp = sumup(a21, a22)
    p2 = square_matrix_multiply_strassens(aTemp, b11)

    bTemp = difference(b12, b22)
    p3 = square_matrix_multiply_strassens(a11, bTemp)

    bTemp = difference(b21, b11)
    p4 = square_matrix_multiply_strassens(a22, bTemp)

    aTemp = sumup(a11, a12)
    p5 = square_matrix_multiply_strassens(aTemp, b22)

    aTemp = difference(a21, a11)
    bTemp = sumup(b11, b12)
    p6 = square_matrix_multiply_strassens(aTemp, bTemp)

    aTemp = difference(a12, a22)
    bTemp = sumup(b21, b22)
    p7 = square_matrix_multiply_strassens(aTemp, bTemp)

    aTemp = sumup(p1, p4)
    bTemp = sumup(aTemp, p7)
    c11 = difference(bTemp, p5)
    c12 = sumup(p3, p5)
    c21 = sumup(p2, p4)

    aTemp = sumup(p1, p3)
    bTemp = sumup(aTemp, p6)
    c22 = difference(bTemp, p2)

    C = [[0 for i in range(0, n)] for j in range(0, n)]
    for i in range(0, new_n):
        for j in range(0, new_n):
            C[i][j] = c11[i][j]
            C[i][j + new_n] = c12[i][j]
            C[i + new_n][j] = c21[i][j]
            C[i + new_n][j + new_n] = c22[i][j]
    return C


print('Strassens Multiplicatioon')
print(square_matrix_multiply_strassens([[0, 1],
                                        [1, 1]], [[0, 1], [1, 1]]))

# ==============================================
# Calculate the power of a matrix in O(k)


def power_of_matrix_navie(A, k):
    """
    Return A^k.
    time complexity = O(k)
    """
    # TODO

    if k == 0:
        return 1
    else:
        s = []
        base_matrix = A
        for i in range(0, k-1):
            if s == []:
                s = square_matrix_multiply_strassens(base_matrix, base_matrix)
            else:
                s = square_matrix_multiply_strassens(s, base_matrix)
        return s


print('Power of matrix_Navie')
print(power_of_matrix_navie([[0, 1], [1, 1]], 3))
# ===================================================
# Calculate the power of a matrix in O(log k)


def power_of_matrix_divide_and_conquer(A, k):
    """
    Return A^k.
    time complexity = O(log k)
    """
    # TODO

    if k == 2:
        return square_matrix_multiply_strassens(A, A)
    elif k == 1:
        return A
    else:
        T1 = power_of_matrix_divide_and_conquer(A, k//2)
        T2 = power_of_matrix_divide_and_conquer(A, k//2 + 1)
        if k % 2 == 0:
            return square_matrix_multiply_strassens(T1, T1)
        else:
            return square_matrix_multiply_strassens(T1, T2)


print('Power of matrix_divide and conquer')
print(power_of_matrix_divide_and_conquer([[0, 1], [1, 1]], 3))

# =================================================


def test():
    assert(find_significant_energy_increase_brute(ENERGY_LEVEL) == (7, 11))
    assert(find_significant_energy_increase_recursive(ENERGY_LEVEL) == (7, 11))
    assert(find_significant_energy_increase_iterative(ENERGY_LEVEL) == (7, 11))
    assert((square_matrix_multiply_strassens([[0, 1], [1, 1]], [
        [0, 1], [1, 1]]) == asarray([
            [1, 1], [1, 2]])).all())
    assert((power_of_matrix_navie([[0, 1], [1, 1]], 3) == asarray([
        [1, 2], [2, 3]])).all())
    assert((power_of_matrix_divide_and_conquer([
        [0, 1], [1, 1]], 3) == asarray([[1, 2], [2, 3]])).all())
    # TODO: Test all of the methods and print results.
    if __name__ == '__main__':
        test()
# ====================================================
