import sys
import time
import itertools


def euclid(a, b):
    r = a % b
    while r != 0:
        a, b = b, r
        r = a % b
    return b


def beuclid(a, b):
    g = 1
    while a & 1 == 0 and b & 1 == 0:
        a >>= 1
        b >>= 1
        g <<= 1

    u, v = a, b
    while u != 0:
        while u & 1 == 0:
            u >>= 1
        while v & 1 == 0:
            v >>= 1
        if u >= v:
            u = u - v
        else:
            v = v - u

    return g * v


def eeuclid(a, b):
    r0, r1 = a, b
    x0, x1 = 1, 0
    y0, y1 = 0, 1

    q, r2 = r0 // r1, r0 % r1
    while r2 != 0:
        x1, x0 = x0 - q * x1, x1
        y1, y0 = y0 - q * y1, y1

        r1, r0 = r2, r1
        q, r2 = r0 // r1, r0 % r1

    return r1, x1, y1


def test_fast():
    # print('Тест быстродействия для пары (78521, 96523)')  # GCD = 1
    print('Тест быстродействия для пары (2852, 4588)')  # GCD = 124

    a, b = 78521, 96523

    start = time.time()
    for _idx in range(1000000):
        euclid(a, b)
    print('Алгоритм Евклида: {}'.format(time.time() - start))

    start = time.time()
    for _idx in range(1000000):
        beuclid(a, b)
    print('Бинарный алгоритм Евклида: {}'.format(time.time() - start))

    start = time.time()
    for _idx in range(1000000):
        eeuclid(a, b)
    print('Расширенный алгоритм Евклида: {}'.format(time.time() - start))


def test_range():
    # print('Тесты в диапазоне (2343, 3343)')
    print('Тесты в диапазоне (10236, 11236)')

    trange = range(2343, 3343)

    start = time.time()
    for (a, b) in itertools.combinations(trange, 2):
        euclid(a, b)
    print('Алгоритм Евклида: {}'.format(time.time() - start))

    start = time.time()
    for (a, b) in itertools.combinations(trange, 2):
        beuclid(a, b)
    print('Бинарный алгоритм Евклида: {}'.format(time.time() - start))

    start = time.time()
    for (a, b) in itertools.combinations(trange, 2):
        eeuclid(a, b)
    print('Расширенный алгоритм Евклида: {}'.format(time.time() - start))


def main():
    assert len(sys.argv) == 4 or (len(sys.argv) == 2 and sys.argv[1] == '-t'), 'Нужен параметр: -e, -b, -ext, -t'
    arg = sys.argv[1]

    if arg == '-t':
        test_fast()
        print('----------------------------')
        test_range()
        sys.exit(0)

    names = {'-e': 'Алгоритм',
             '-b': 'Бинарный алгоритм',
             '-ext': 'Расширенный алгоритм'}

    print('{} Евклида'.format(names[arg]))
    a = int(sys.argv[2])
    b = int(sys.argv[3])

    assert a > 0 and b > 0, 'Числа должны быть положительными'
    a, b = max(a, b), min(a, b)

    print('-------------------------')
    if arg == '-e':
        print('НОД({}, {}) = {}'.format(a, b, euclid(a, b)))
    elif arg == '-b':
        print('НОД({}, {}) = {}'.format(a, b, beuclid(a, b)))
    elif arg == '-ext':
        d, x, y = eeuclid(a, b)
        print('НОД({}, {}) = {}'.format(a, b, d))
        print('{} * ({}) + {} * ({}) = {}'.format(a, x, b, y, d))
    else:
        raise Exception('Неверное значение параметра: -e, -b, -ext')


if __name__ == '__main__':
    main()
