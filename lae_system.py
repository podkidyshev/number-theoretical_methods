from utils import Vector as V, Matrix as M, ratio


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
    while True:
        if not V.is_zero(ws[-1]) and mul_a(ws[-1], ws[-1], les.a, p) == 0:
            return
        if V.is_zero(ws[-1]):
            x = V.zero(les.n)
            for wi in ws[:-1]:
                scalar = ratio(V.mul_sum(wi, ws[0], p),
                               V.mul_sum(wi, wi, p), p)
                x = V.add(x, V.mul_scalar(wi, scalar, les.p), p)
            return x

        vec = V.zero(les.n)
        for wj in ws:
            alpha = alpha_ij(ws[-1], wj, les.a, p)
            vec = V.add(vec, V.mul_scalar(wj, alpha, les.p), p)
        vec = V.mul_scalar(vec, -1, p)
        wi = V.minus(M.mul_vec(les.a, ws[-1], p), vec, p)
        ws.append(wi)


def lanczos2(les: LinearEquationSystem):
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


if __name__ == '__main__':
    test_system = [
        '1 1 3 5',
        '1 6 4 4',
        '3 4 5 5']  # 1, 2, 3 is solution too
    test_p = 7
    test_les = LinearEquationSystem(test_system, test_p)
    res = lanczos2(test_les)
    print(res)
    print(M.mul(test_les.a, M.t([res]), test_les.p) == test_les.b)
