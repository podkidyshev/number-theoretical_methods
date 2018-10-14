import sys
import math
import random

import utils
import prime
import factor
import euclid


def factorize(n):
    base = {}
    n, k = utils.fac2k(n)
    if k > 0:
        base[2] = k
    while n > 1:
        divisor = factor.ppollard(n) if not prime.isprime(n) else n
        if divisor is None:
            raise ValueError('Ошибка при факторизации')
        if divisor not in base:
            base[divisor] = 0
        base[divisor] += 1
        n //= divisor
    return base


def primitive_root(p, max_count=0):
    """Только для простых p"""
    assert prime.isprime(p)
    degs = [(p - 1) // f for f in factorize(p - 1)]  # факторизация фи(p) = p-1
    answer = []
    for g in range(1, p):
        print(g, end='\r')
        if all([pow(g, d, p) != 1 for d in degs]):
            if not max_count:
                return g
            answer.append(g)
    return answer[:max_count]


def square_root(n):
    x1 = n
    x2 = int((x1 + (n / x1)) / 2)
    while x2 < x1:
        x1, x2 = x2, int((x2 + (n / x2)) / 2)
    return x1


def shanks(m, g, h):
    r = square_root(m) + 1
    pairs = {pow(g, a, m): a for a in range(r)}
    g1 = pow(utils.get_inverse(g, m), r, m)
    for b in range(r):
        value = (pow(g1, b, m) * h) % m
        if value in pairs:
            return pairs[value] + r * b


def equation(a, b, m):
    # ax = b (mod m)
    a, b = a % m, b % m
    d = euclid.euclid(a, m)  # число решений
    if b % d != 0:
        return  # нет решений, если d не делит b

    a_new, b_new, m_new = a // d, b // d, m // d
    d_new, q, r = euclid.eeuclid(a_new, m_new)  # 1 = a * q + m * r
    q, r = q % m, r % m

    x0 = (b_new * q) % m_new
    for j in range(d):
        yield x0 + m_new * j


def next_yab(y, a, b, params):
    g, h, m = params
    if y <= m // 3:
        y = (y * h) % m
        a = (a + 1) % (m - 1)
    elif m // 3 < y <= 2 * m // 3:
        y = (y * y) % m
        a = (a * 2) % (m - 1)
        b = (b * 2) % (m - 1)
    elif 2 * m // 3 < m:
        y = (y * g) % m
        b = (b + 1) % (m - 1)
    return y, a, b


def ppollard(m, g, h, e=0.05):
    """http://math.ucdenver.edu/graduate/thesis/dreiff.pdf"""
    t = square_root(2 * m * math.log(1 / e)) + 1
    while True:
        i = 1
        s = random.randint(0, m - 2)

        yi, ai, bi = next_yab(pow(g, s, m), 0, s, (g, h, m))
        y2i, a2i, b2i = next_yab(yi, ai, bi, (g, h, m))

        while i < t and yi != y2i:
            i += 1
            yi, ai, bi = next_yab(yi, ai, bi, (g, h, m))
            y2i, a2i, b2i = next_yab(*next_yab(y2i, a2i, b2i, (g, h, m)), (g, h, m))

        if yi == y2i:
            aa, bb = (a2i - ai) % (m - 1), (bi - b2i) % (m - 1)  # h^aa = g^bb
            d = euclid.euclid(aa, m - 1)
            if d < square_root(m - 1):
                for x in equation(aa, bb, m - 1):
                    if pow(g, x, m) == h:
                        return x


FUNCS = [(shanks, 'Алгоритм Шенкса'), (ppollard, 'p-метод Полларда')]


def test_accuracy(m, g, h):
    print('Тест корректности алгоритмов по нахождению решения g^x = h\ng - образующая циклической группы порядка m')
    print('Образующий g =', g)
    print('Порядок подгруппы m =', m)
    print('Значение степени h =', h)

    for (func, func_name) in FUNCS:
        x = func(m, g, h)
        assert h == pow(g, x, m)
        print('{}: {}'.format(func_name, x if x is not None else 'ответ не найден'))


if __name__ == '__main__':
    if sys.argv[1] == '-a':
        test_m = int(sys.argv[2]) if len(sys.argv) > 2 else 29
        test_g = int(sys.argv[3]) % test_m if len(sys.argv) > 2 else 2
        test_h = int(sys.argv[4]) % test_m if len(sys.argv) > 2 else 17
        test_accuracy(test_m, test_g, test_h)
    elif sys.argv[1] == '-p':
        test_m = int(sys.argv[2]) if len(sys.argv) > 2 else 29
        if '-all' in sys.argv:
            print('Все образующие группы порядка m = {} (выведем макимум 100 первых):'.format(test_m))
            print(primitive_root(test_m, 100))
        else:
            print('Первый образующий элемент группы порядка m = {}:'.format(test_m), primitive_root(test_m))

# 105251
# 55609
