import sys
import random

from utils import Vector as V, Matrix as M, Polynomial as Poly, ratio


class LinearEquationSystem:
    def __init__(self, matrix, p):
        matrix = [row.split() for row in matrix]
        self.p = p                   # поле вычетов
        self.n = len(matrix)         # число уравнений
        self.m = len(matrix[0]) - 1  # число неизвестных

        self.a = [[int(matrix[i][j]) for j in range(self.m)] for i in range(self.n)]
        self.b = [[int(matrix[i][-1])] for i in range(self.n)]


def gaussian(les: LinearEquationSystem):
    pass


def mul_a(x, y, a, p):
    return V.mul_sum(M.mul_vec(a, x, p), y, p)


def alpha_ij(wi, wj, a, p):
    nom = mul_a(M.mul_vec(a, wi, p), wj, a, p)
    rat = mul_a(wj, wj, a, p)
    return ratio(nom, rat, p)


def lanczos(les: LinearEquationSystem):
    assert les.n == les.m
    p = les.p
    ws = [M.column(les.b, 0)]
    vs = [None, M.mul_vec(les.a, ws[0], p)]
    ws.append(V.minus(vs[1],
                      V.mul_scalar(ws[0],
                                   ratio(V.mul_sum(vs[1], vs[1], p),
                                         V.mul_sum(ws[0], vs[1], p), p), p), p))

    while True:
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
    while True:
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
def gen_test(size, p):
    a = [[random.randint(0, p - 1) for _j in range(size)] for _i in range(size)]
    while M.det(a, p) == 0:
        a = [[random.randint(0, p - 1) for _j in range(size)] for _i in range(size)]
    b = [[random.randint(0, p - 1)] for _i in range(size)]
    for row, bi in zip(a, b):
        print(*row, *bi)


if __name__ == '__main__':
    operation = sys.argv[1]

    if operation == '-a':
        # _system = [
        #     '1 1 3 5',
        #     '1 6 4 4',
        #     '3 4 6 5'
        # ]
        # _system = [
        #     '1 0 4 0 3',
        #     '0 6 0 1 4',
        #     '0 1 2 0 3',
        #     '0 1 0 3 1'
        # ]
        # _system = [
        #     '0 0 1 1 2 4 0 3',
        #     '4 4 1 6 6 0 6 0',
        #     '2 6 3 6 6 4 1 3',
        #     '2 2 6 5 3 3 1 0',
        #     '0 6 5 5 6 6 6 4',
        #     '4 0 3 2 5 2 5 2',
        #     '6 6 5 6 0 2 4 4'
        # ]
        _system = [
            '12 17 9 5 13 24 22 2 9 5 15',
            '4 19 5 21 2 26 20 21 8 21 17',
            '16 21 1 28 21 24 16 16 7 8 15',
            '24 25 12 26 9 27 27 17 5 19 15',
            '23 15 15 18 0 1 19 3 2 16 10',
            '12 16 5 11 27 11 4 15 1 2 11',
            '0 24 14 6 26 6 25 8 7 17 21',
            '20 6 24 20 17 6 0 25 7 13 9',
            '5 7 14 24 9 28 15 8 11 8 0',
            '5 7 13 18 6 17 5 8 24 10 0'
        ]
        _p = 29  # 349
        _les = LinearEquationSystem(_system, _p)
        assert M.det(_les.a, _les.p) != 0
        _res = wiedemann2(_les)
        print(_res)
        print(M.mul(_les.a, M.t([_res]), _les.p) == _les.b)
    elif operation == '-g':
        gen_test(10, 29)
    else:
        print('Wrong operation')
