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
        q, r = Poly.ratio(r0, r1, p)
        v = Poly.minus(v0, Poly.mul(q, v1, p), p)
        v0, v1 = v1, v
        r0, r1 = r1, r

    v1, r2 = Poly.ratio(v1, [v1[-1]], p)
    assert r2 == [0]
    assert Poly.compute2(v1, seq[:len(v1)], p) == 0
    assert len(v1) <= n + 1
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


def wiedemann(les: LinearEquationSystem):
    a, b, p, n = les.a, M.t(les.b)[0], les.p, les.n
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

        assert len(seq)
        f = berlekamp(seq, p)

        ys.append(V.add(ys[k], M.mul_vec(f_tilda(f, a, p), b, p), p))
        bs.append(V.add(b, M.mul_vec(a, ys[k + 1], p), p))
        ds.append(ds[k] + len(f))
        k += 1
    return V.mul_scalar(ys[k], -1, p)


if __name__ == '__main__':
    _system = [
        '1 1 3 5',
        '1 6 4 4',
        '3 4 6 5']
    _p = 7
    _les = LinearEquationSystem(_system, _p)
    _res = wiedemann(_les)
    print(_res)
    print(M.mul(_les.a, M.t([_res]), _les.p) == _les.b)
    # print(berlekamp([1, 6, 1, 6, 1, 6], 7))
    # print(berlekamp([2, 1, 1, 5, 6, 6], 7))  # 3x + 2y
    # print(berlekamp([0, 1, 2, 4, 6, 0], 7))
