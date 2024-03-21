from numpy import asarray

def sum_up(x,y):
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

    assert (len(x) & (len(x) -1)) == 0, "A is not a power of 2"
    n = len(x)

    if n == 1:
        C = [[0 for j in range(0, n)] for i in range(0, n)]
        for i in range(0, n):
            for j in range(0, n):
                C[i][j] = x[i][j] * y[i][j]
        return C
    else:  # dividing the input matrices A and B
        new_n = int(n / 2)

    a11 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]
    a12 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]
    a21 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]
    a22 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]

    b11 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]
    b12 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]
    b21 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]
    b22 = [[0 for i in range(0, new_n)] for j in range(0, new_n)]

    aTemp = [[0 for i in range(0, new_n)] for j in range(0, new_n)]
    bTemp = [[0 for i in range(0, new_n)] for j in range(0, new_n)]

    for i in range(0, new_n):
        for j in range(0, new_n):
            a11[i][j] = x[i][j]
            a12[i][j] = x[i][j + new_n]
            a21[i][j] = x[i + new_n][j]
            a22[i][j] = x[i + new_n][j + new_n]

            b11[i][j] = y[i][j]
            b12[i][j] = y[i][j + new_n]
            b21[i][j] = y[i + new_n][j]
            b22[i][j] = y[i + new_n][j + new_n]

    aTemp = sum_up(a11, a22)
    bTemp = sum_up(b11, b22)
    p1 = square_matrix_multiply_strassens(aTemp, bTemp)

    aTemp = sum_up(a21, a22)
    p2 = square_matrix_multiply_strassens(aTemp, b11)

    bTemp = difference(b12, b22)
    p3 = square_matrix_multiply_strassens(a11, bTemp)

    bTemp = difference(b21, b11)
    p4 = square_matrix_multiply_strassens(a22, bTemp)

    aTemp = sum_up(a11, a12)
    p5 = square_matrix_multiply_strassens(aTemp, b22)

    aTemp = difference(a21, a11)
    bTemp = sum_up(b11, b12)
    p6 = square_matrix_multiply_strassens(aTemp, bTemp)

    aTemp = difference(a12, a22)
    bTemp = sum_up(b21, b22)
    p7 = square_matrix_multiply_strassens(aTemp, bTemp)

    aTemp = sum_up(p1, p4)
    bTemp = sum_up(aTemp, p7)
    c11 = difference(bTemp, p5)
    c12 = sum_up(p3, p5)
    c21 = sum_up(p2, p4)

    aTemp = sum_up(p1, p3)
    bTemp = sum_up(aTemp, p6)
    c22 = difference(bTemp, p2)

    C = [[0 for i in range(0, n)] for j in range(0, n)]
    for i in range(0, new_n):
        for j in range(0, new_n):
            C[i][j] = c11[i][j]
            C[i][j + new_n] = c12[i][j]
            C[i + new_n][j] = c21[i][j]
            C[i + new_n][j + new_n] = c22[i][j]
    return C


########################################################################################

def power_of_matrix_navie(A, k):
    if k==0:
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


########################################################################################

def power_of_matrix_divide_and_conquer(A, k):
    if k==2:
        return square_matrix_multiply_strassens(A, A)
    elif k == 1:
        return A
    else:
        T1 = power_of_matrix_divide_and_conquer(A, k//2)
        T2 = power_of_matrix_divide_and_conquer(A, k//2 + 1 )
        if k%2 == 0:
            return square_matrix_multiply_strassens(T1, T1)
        else:
            return square_matrix_multiply_strassens(T1, T2)


########################################################################################


assert((square_matrix_multiply_strassens([[0, 1], [1, 1]], [[0, 1], [1, 1]]) ==
        asarray([[1, 1], [1, 2]])).all())
assert((power_of_matrix_navie([[0, 1], [1, 1]], 3) ==
                                    asarray([[1, 2], [2, 3]])).all())
assert((power_of_matrix_divide_and_conquer([[0, 1], [1, 1]], 3) ==
                                               asarray([[1, 2], [2, 3]])).all())

print(square_matrix_multiply_strassens([[0, 1], [1, 1]], [[0, 1], [1, 1]]))
print(power_of_matrix_navie([[0, 1], [1, 1]], 3))
print(power_of_matrix_divide_and_conquer([[0, 1], [1, 1]], 3))