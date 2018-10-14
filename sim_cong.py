# simultaneous congruences
import sys
import functools
import itertools

import euclid


def inverse(a, m):  # in ring
    a %= m
    d, x, y = euclid.eeuclid(a, m)  # a * x + m * y = 1 (mod m)
    if d != 1:
        raise Exception('Нет обратного для {} и {}'.format(a, m))
    return x % m


def chinese_remainder_theorem(congs):
    m = functools.reduce(int.__mul__, (cong[1] for cong in congs))
    u = 0
    for cong in congs:
        c = m // cong[1]
        d = inverse(m // cong[1], cong[1])
        u = (u + cong[0] * c * d) % m
    return u


def garner(congs):
    c = [1]
    for i in range(1, len(congs)):
        c.append(1)
        for j in range(0, i):
            u = inverse(congs[j][1], congs[i][1])
            c[i] = (u * c[i]) % congs[i][1]

    u = congs[0][0]
    x = u

    for i in range(1, len(congs)):
        u = ((congs[i][0] - x) * c[i]) % congs[i][1]
        x += u * functools.reduce(int.__mul__, [congs[j][1] for j in range(0, i)])
    return x


def main():
    assert len(sys.argv) == 2, 'Нужен параметр: -gc, -ga'
    arg = sys.argv[1]

    names = {'-gc': ' с помощью греко-китайской теоремы',
             '-ga': 'методом Гарнера'}

    print('Решить систему сравнений {}'.format(names[arg]))
    n = int(input('Введите число уравнений в системе: '))
    print('Введите {} пар чисел a_k и m_k через пробел на отдельных строках при x ≡ a_k (mod m_k): '.format(n))
    congs = [list(map(int, input().split())) for _idx in range(n)]

    for comb in itertools.combinations((cong[1] for cong in congs), 2):
        if euclid.euclid(*comb) != 1:
            raise Exception('Пара не взаимно простых m_k: {} и {}'.format(comb[0][1], comb[1][1]))

    for cong in congs:
        cong[0] %= cong[1]

    if arg == '-gc':
        print('Решение: x = {}'.format(chinese_remainder_theorem(congs)))
    elif arg == '-ga':
        print('Решение: x = {}'.format(garner(congs)))
    else:
        raise Exception('Неверное значение параметра: -gc, -ga')


if __name__ == '__main__':
    main()

"""
Тесты:
-----------------------------------------------------------------
http://pmpu.ru/vf4/modular/crt
3
7 8
-1 11
3 15

x = 1143
-----------------------------------------------------------------
https://ru.wikipedia.org/wiki/Китайская_теорема_об_остатках
3
1 2
2 3
6 7

x = 41
-----------------------------------------------------------------
"""
