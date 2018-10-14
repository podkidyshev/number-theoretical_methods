import sys
import time
from math import log, sqrt, exp
import random
import functools
import sympy

import utils
import euclid
import gaussian
from chain_fractions import *

B_DEFAULT = (2, 3, 5)
N_DEFAULT = 1728239


def compressor(x):
    return x * x + 1


def ppollard(n, c=1, f=compressor):
    a = c
    b = c
    while True:
        a, b = f(a) % n, f(f(b) % n) % n
        d = euclid.euclid(a - b, n)
        if 1 < d < n:
            return d
        elif d == n:
            return


def p1pollard(n):
    b = p1pollard_base(n)
    a = random.randint(2, n - 2)
    d = euclid.euclid(a, n)
    if d >= 2:
        return d
    for pi in b:
        li = int(log(n) / log(pi))
        a = pow(a, pow(pi, li), n)
    d = euclid.euclid(a - 1, n)
    if 1 < d < n:
        return d


def p1pollard_base(n):
    base_size = int(n ** (1 / 6))
    return utils.generate_base(min(base_size, utils.MAX_B_SIZE))


def is_b_smooth(p, b):
    alpha = []
    for bi in b:
        k = 0
        while p % bi == 0 and not bi < 0 < p and p != 1:
            p //= bi
            k += 1
        alpha.append(k)
    return p == 1, alpha, [al % 2 for al in alpha]


def dixon(n):
    base = dixon_base(n)
    base = [-1] + list(filter(lambda bi: utils.legendre(n, bi) == 1, base))
    h = len(base) - 1

    ps = []
    alphas = []
    es = []
    convergent = gen_convergent(gen_square_chain_fraction(n))
    while len(ps) < h + 2:
        try:
            pi, qi = next(convergent)
        except ValueError:
            return None

        pi2 = pi ** 2 % n
        if n - pi2 < pi2:
            pi2 = -(n - pi2)
        smooth, alpha, e = is_b_smooth(pi2, base)
        if smooth:
            ps.append(pi)
            alphas.append(alpha)
            es.append(e)

    for ks in gaussian.gen_gaussian(es):
        s = 1
        for k in ks:
            s = (s * ps[k]) % n
        t = 1
        for b_idx, b in enumerate(base):
            t = (t * pow(b, functools.reduce(int.__add__, (alphas[k][b_idx] for k in ks)) // 2, n)) % n
        # проверка, что ks - не решение системы
        assert pow(s, 2, n) == pow(t, 2, n)

        if s != t and s != n - t:
            p = euclid.euclid((s - t) % n, n)
            return p


def dixon_base(n):
    base_size = int(sqrt(exp(2 / 3 * sqrt(log(n) * log(log(n))))))
    return utils.generate_base(min(base_size, utils.MAX_B_SIZE))


FUNCS = [('p-метод Полларда', ppollard),
         ('(p-1)-метод Полларда', p1pollard),
         ('Метод непрерывных дробей', dixon)]


def test_accuracy(n):
    print('Тест корректности числа {}'.format(n))
    for (func_name, func) in FUNCS:
        divisor = func(n)
        print(func_name + ':', 'Делитель не найден' if divisor is None else divisor)


def test_speed(n):
    print('Тест скорости работы алгоритмов (10000 запусков) на числе {}:'.format(n))
    divisors = [str(func(n)) for (func_name, func) in FUNCS]
    times = []
    for (func_name, func) in FUNCS:
        start = time.time()
        for _idx in range(10000):
            func(n)
        times.append(time.time() - start)

    print(' ' * 10, 'p-метод Полларда  ', '(p-1)-метод Полларда  ', 'Метод непрерывных дробей')
    print('{:>9}   {:>15}   {:>20}   {:>23}'.format('Делитель', *divisors))
    print('{:>9}   {:>13.3f} c   {:>18.3f} c   {:>21.3f} c'.format('Время', *times))


if __name__ == '__main__':
    arg = sys.argv[1]
    test_n, k = utils.fac2k(int(sys.argv[2]))
    if sympy.isprime(test_n):
        print('{} - простое'.format(test_n))
        exit(0)
    if arg == '-a':
        test_accuracy(test_n)
    elif arg == '-s':
        test_speed(test_n)

# 141642509 * 179425457 = 25414271907951613
# 63169 - составное
# 21299881 - составное
# 37089829131812000665041223352794276544031149 - очень составное
