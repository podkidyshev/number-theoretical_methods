from functools import reduce
from operator import add

import utils


class LinearEquationSystem:
    def __init__(self, matrix, p):
        matrix = [row.split() for row in matrix]
        self.p = p                   # поле вычетов
        self.n = len(matrix)         # число уравнений
        self.m = len(matrix[0]) - 1  # число неизвестных

        self.a = [[int(matrix[i][j]) for j in range(len(matrix[0]) - 1)] for i in range(len(matrix))]
        self.b = [[int(matrix[i][-1])] for i in range(len(matrix))]


def gaussian(les: LinearEquationSystem):
    pass


def vector_zero(x):
    return all(el == 0 for el in x)


def vector_sum(x, y, p):
    return [(xx + yy) % p for xx in x for yy in y]


def vector_mul(x, alpha, p):
    return [(xx * alpha) % p for xx in x]


def matrix_transpose(a):
    return [[a[i][j] for i in range(len(a))] for j in range(len(a[0]))]


def matrix_column(a, j):
    return [a[i][j] for i in range(len(a))]


def matrix_mul(a, b, p):
    assert len(a[0]) == len(b)
    return [[reduce(add, scalar_mul(a[i], matrix_column(b, j), p)) % p for j in range(len(b[0]))] for i in range(len(a))]


def scalar_mul(x, y, p):
    assert len(x) == len(y)
    return [(x[i] * y[i]) % p for i in range(len(x))]


def scalar_mul_a(x, y, a, p):
    x = matrix_transpose([x])
    ax = matrix_transpose(matrix_mul(a, x, p))[0]
    return reduce(add, scalar_mul(ax, y, p)) % p


def mul_matvec(a, x, p):
    # res = Ax
    return matrix_column(matrix_mul(a, matrix_transpose([x]), p), 0)


def alpha_ij(wi, wj, les: LinearEquationSystem):
    nom = scalar_mul_a(mul_matvec(les.a, wi, les.p), wj, les.a, les.p)
    rat = scalar_mul_a(wj, wj, les.a, les.p)
    return (nom * utils.get_inverse(rat, les.p)) % les.p


def lanczos(les: LinearEquationSystem):
    assert les.n == les.m
    ws = [matrix_column(les.b, 0)]
    while True:
        if not vector_zero(ws[-1]) and scalar_mul_a(ws[-1], ws[-1], les.a, les.p) == 0:
            return
        if vector_zero(ws[-1]):
            pass  # SOLUTION
        vec = [0] * les.n
        for j in range(len(ws)):
            vec = vector_sum(vec, vector_mul(ws[j], alpha_ij(ws[-1], ws[j], les), les.p), les.p)
        wi = vector_sum(mul_matvec(les.a, ws[-1], les.p), vector_mul(vec, -1, les.p), les.p)
        ws.append(wi)


if __name__ == '__main__':
    # test_a = [[1, 1, 1],
    #           [1, 1, 1],
    #           [1, 1, 1]]
    test_b = [[1, 2, 1],
              [2, 2, 2],
              [1, 2, 0]]
    test_x = [1, 1, 1]
    test_y = [2, 3, 4]
    res = scalar_mul_a(test_x, test_y, test_b, 3)
    # utils.matrix_print([res])
    print(res)
