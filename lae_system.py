import copy
import sys
import random

from utils import Vector as V, Matrix as M, Polynomial as Poly, ratio, get_inverse

MAX_TRIALS = 20


class LinearEquationSystem:
    def __init__(self, matrix, p):
        matrix = [row.split() for row in matrix]
        self.p = p  # поле вычетов
        self.n = len(matrix)  # число уравнений
        self.m = len(matrix[0]) - 1  # число неизвестных

        self.a = [[int(matrix[i][j]) for j in range(self.m)] for i in range(self.n)]
        self.b = [[int(matrix[i][-1])] for i in range(self.n)]


def find_nonzero(a, col):
    for row in range(col, len(a)):
        if a[row][col] != 0:
            return row
    return -1


def gaussian(les: LinearEquationSystem):
    a, b = copy.deepcopy(les.a), copy.deepcopy(les.b)
    n, m, p = les.n, les.m, les.p

    for col in range(m):
        row_replaced = find_nonzero(a, col)
        if row_replaced == -1:
            continue
        if row_replaced != col:
            a[col], a[row_replaced] = a[row_replaced], a[col]
            b[col], b[row_replaced] = b[row_replaced], b[col]
        for row in range(col + 1, n):
            if a[row][col] != 0:
                scalar = get_inverse(a[row][col], p) * a[col][col] % p
                a[row] = V.mul_scalar(a[row], scalar, p)
                a[row] = V.minus(a[row], a[col], p)
                brow = b[row][0]
                b[row] = [(brow * scalar - b[col][0]) % p]
                assert a[row][col] == 0

    x = [0] * m
    for row in range(m - 1, -1, -1):
        xcol = b[row][0]
        for col in range(m - 1, row, -1):
            xcol = (xcol - a[row][col] * x[col]) % p
        xcol = xcol * get_inverse(a[row][row], p) % p
        x[row] = xcol
    return x


def mul_a(x, y, a, p):
    return V.mul_sum(M.mul_vec(a, x, p), y, p)


def alpha_ij(wi, wj, a, p):
    nom = mul_a(M.mul_vec(a, wi, p), wj, a, p)
    rat = mul_a(wj, wj, a, p)
    return ratio(nom, rat, p)


def lanczos(les: LinearEquationSystem):
    assert les.n == les.m
    for row in range(les.n):
        for col in range(row + 1):
            if les.a[row][col] != les.a[col][row]:
                print('Алгоритм Ланцоша : определен только для симмтеричных матриц, '
                      'см. индексы ({0},{1}) и ({1},{0})'.format(row + 1, col + 1))
                return
    p = les.p
    ws = [M.column(les.b, 0)]
    vs = [None, M.mul_vec(les.a, ws[0], p)]
    ws.append(V.minus(vs[1],
                      V.mul_scalar(ws[0],
                                   ratio(V.mul_sum(vs[1], vs[1], p),
                                         V.mul_sum(ws[0], vs[1], p), p), p), p))

    for trial in range(max(MAX_TRIALS, les.n)):
        if V.mul_sum(ws[-1], M.mul_vec(les.a, ws[-1], p), p) == 0:
            if V.is_zero(ws[-1]):
                x = V.zero(les.n)
                for i in range(len(ws) - 1):
                    bi = ratio(V.mul_sum(ws[i], ws[0], p),
                               V.mul_sum(ws[i], vs[i + 1], p), p)
                    x = V.add(x, V.mul_scalar(ws[i], bi, p), p)
                return x

        vs.append(M.mul_vec(les.a, ws[-1], p))
        w1 = V.mul_scalar(ws[-1], ratio(V.mul_sum(vs[-1], vs[-1], p),
                                        V.mul_sum(ws[-1], vs[-1], p), p), p)
        w2 = V.mul_scalar(ws[-2], ratio(V.mul_sum(vs[-1], vs[-2], p),
                                        V.mul_sum(ws[-2], vs[-2], p), p), p)
        w_end = V.minus(V.minus(vs[-1], w1, p),
                        w2, p)
        ws.append(w_end)


def berlekamp(seq: list, p: int):
    n = len(seq) // 2
    r0 = [1] + [0] * (2 * n)
    r1 = Poly.shrink(seq[:])
    v0, v1 = [0], [1]

    while n <= len(r1):
        if r1 == [0]:
            return None
        q, r = Poly.ratio(r0, r1, p)
        v = Poly.minus(v0, Poly.mul(q, v1, p), p)
        v0, v1 = v1, v
        r0, r1 = r1, r

    if v1[-1] != 0:
        v1, r2 = Poly.ratio(v1, [v1[-1]], p)
        assert r2 == [0]
    assert all(Poly.compute2(v1, seq[i:i + len(v1)], p) == 0 for i in range(len(seq) - len(v1)))
    return v1


def fa(f, a, p):
    res = M.zero(len(a), len(a))
    for power, fi in enumerate(reversed(f)):
        res = M.sum(res, M.mul_scalar(M.power(a, power, p), fi, p), p)
    return res


def f_tilda(f, a, p):
    gz = fa(f, a, p)
    g0 = fa(f, M.zero(len(a), len(a)), p)
    numenator = M.sum(gz, M.mul_scalar(g0, p - 1, p), p)
    return M.mul(numenator, M.inverse(a, p), p)


def wiedemann1(les: LinearEquationSystem):
    a, b, p, n = les.a, M.t(les.b)[0], les.p, les.n
    for trial in range(MAX_TRIALS):
        bs = [b]
        ys = [[0] * n]
        ds = [0]
        k = 0

        while bs[k] != V.zero(n):
            u = [random.randint(0, p - 1) for _idx in range(n)]
            seq = []
            for i in range(2 * (n - ds[k])):
                ai = M.power(a, i, p)
                aib = M.mul_vec(ai, b, p)
                uaib = V.mul_sum(u, aib, p)
                seq.append(uaib)

            if not len(seq):
                break
            f = berlekamp(seq, p)
            if f is None:
                break

            ys.append(V.add(ys[k], M.mul_vec(f_tilda(f, a, p), b, p), p))
            bs.append(V.add(b, M.mul_vec(a, ys[k + 1], p), p))
            ds.append(ds[k] + len(f))
            k += 1
        if bs[k] != V.zero(n):
            print(trial, end='\r')
            continue
        return V.mul_scalar(ys[k], -1, p)


def wiedemann2(les: LinearEquationSystem):
    a, b, p, n = les.a, M.t(les.b)[0], les.p, les.n
    a_powers = [M.unit(n, n)]
    for i in range(1, 2 * n):
        a_powers.append(M.mul(a_powers[-1], a, p))
    aib = [M.mul_vec(ai, b, p) for ai in a_powers]
    k = 0
    gs = [[1]]
    uk1 = [0, 1] + ([0] * (n - 2))
    while Poly.deg(gs[k]) < n and k < n:
        seq = []
        for i in range(2 * n - Poly.deg(gs[k])):
            gab = M.mul_vec(fa(gs[k], a, p), aib[i], p)
            ugab = V.mul_sum(uk1, gab, p)
            seq.append(ugab)

        assert len(seq)
        f = berlekamp(seq, p)
        gs.append(Poly.mul(f, gs[k], p))
        k += 1
        uk1 = ([0] * k) + [1] + ([0] * (n - k - 1))

    print('k =', k)
    g = gs[-1]
    x = V.zero(n)
    for i in range(Poly.deg(g)):
        x = V.add(x, V.mul_scalar(aib[i], -g[i], p), p)
    return x


#######################################################################
def gen_test(size, p, sym):
    a = [[random.randint(0, p - 1) for _j in range(size)] for _i in range(size)]
    while M.det(a, p) == 0:
        a = [[random.randint(0, p - 1) for _j in range(size)] for _i in range(size)]
    if sym:
        for i in range(len(a)):
            for j in range(i + 1):
                a[i][j] = a[j][i]
    b = [[random.randint(0, p - 1)] for _i in range(size)]

    print(size, p)
    for row, bi in zip(a, b):
        print(*row, *bi)


def read(filename):
    with open(filename) as f:
        n, p = list(map(int, f.readline().split(' ')))
        a = [f.readline() for _idx in range(n)]
        return p, n, a


#######################################################################
FUNCS = [(gaussian, 'Метод Гаусса'),
         (lanczos, 'Алгоритм Ланцоша'),
         (wiedemann1, 'Алгоритм Видемана (вероятностный)')]


def test_accuracy(les):
    for func, func_name in FUNCS:
        x = func(les)
        if x is None:
            print(func_name, ': ответ не найден')
        else:
            correct = M.mul(les.a, M.t([x]), les.p) == les.b
            print(func_name, ':', x, 'Ответ верный' if correct else 'Ответ неверный')


def main():
    operation = sys.argv[1]

    if operation == '-a':
        _p, _n, _a = read('test.txt')
        _les = LinearEquationSystem(_a, _p)
        test_accuracy(_les)
    elif operation == '-g':
        gen_test(15, 31, '-s' in sys.argv)
    else:
        print('Некорректная операция')


if __name__ == '__main__':
    main()
