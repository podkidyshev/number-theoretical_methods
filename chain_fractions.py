import math


def gen_chain_fraction(p, q):
    a = int(p / q)
    yield a

    while p != q:
        p, q = q, p - q * a
        a = int(p / q)
        yield a


def gen_square_chain_fraction(n):
    a0 = math.sqrt(n)
    r0 = int(a0)
    yield r0
    ratio0 = 1
    numenator0 = 0

    while True:
        numenator1 = r0 * ratio0 - numenator0
        ratio1 = (n - numenator1 * numenator1) // ratio0
        if ratio1 == 0:
            raise ValueError('Число является полным квадратом')
        r1 = int((a0 + numenator1) / ratio1)

        yield r1
        r0, ratio0, numenator0 = r1, ratio1, numenator1


def gen_convergent(generator):
    p0, p1 = 0, 1
    q0, q1 = 1, 0

    while True:
        ai = next(generator)
        pi = ai * p1 + p0
        qi = ai * q1 + q0
        yield pi, qi
        p0, p1, q0, q1 = p1, pi, q1, qi


if __name__ == '__main__':
    pass
