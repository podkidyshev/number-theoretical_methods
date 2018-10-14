import sys
import random
import functools
import time

import utils


def fermat(n):
    a = random.randint(2, n - 2)
    r = pow(a, n - 1, n)
    return r == 1


def solovay_strassen(n):
    a = random.randint(2, n - 2)
    r = pow(a, (n - 1) // 2, n)

    if r != 1 and r != n - 1:
        return False

    s = utils.jacobi(a, n) % n
    return r == s


def miller_rabin(n):
    r, s = utils.fac2k(n - 1)
    a = random.randint(2, n - 2)
    y = pow(a, r, n)

    if y != 1 and y != n - 1:
        j = 1
        while j <= s - 1 and y != n - 1:
            y = (y * y) % n
            if y == 1:
                return False
            j += 1
        if y != n - 1:
            return False
    return True


def isprime(n, rounds=10):
    if n > 5:
        for round in range(rounds):
            if not miller_rabin(n):
                return False
        return True
    return n in [2, 3, 5]


FUNCS = [(fermat, 'Тест Ферма'),
         (solovay_strassen, 'Тест Соловея-Штрассена'),
         (miller_rabin, 'Тест Миллера-Рабина')]
RES_TABLE = {True: 'простое',
             False: 'составное'}


def test_accuracy(n, rounds):
    print('Тестируем число {} на {} раундах'.format(n, rounds))
    print('-' * 73)
    for (test_func, test_func_name) in FUNCS:
        res = all(test_func(n) for _idx in range(rounds))
        print('{}: {}'.format(test_func_name, RES_TABLE[res]))


def test_statistics(n):
    print('Тестируем количество ошибок для числа {} на 1000 запусков программы'.format(n))

    table = []
    for rounds in range(1, 5):
        table_rounds = []
        for (test_func, test_func_name) in FUNCS:
            errors = functools.reduce(int.__add__, [int(all(test_func(n) for _round in range(rounds))) for _idx in range(1000)])
            table_rounds.append(errors)
        table.append(table_rounds)

    print('Число раундов  ', 'Тест Ферма  ', 'Тест Соловея-Штрассена  ', 'Тест Миллера-Рабина')
    for rounds, table_rounds in enumerate(table):
        print('{:>13}   {:>10}   {:>22}   {:>19}'.format(rounds + 1, *tuple(table_rounds)))


def test_speed(n):
    print('Тестируем быстродействие алгоритмов на числе {}'.format(n))

    for (test_func, test_func_name) in FUNCS:
        start = time.time()
        for _idx in range(1000000):
            test_func(n)
        start = time.time() - start
        print('{}: {:.2f} c'.format(test_func_name, start))


if __name__ == '__main__':
    if sys.argv[1] == '-a':
        test_number = int(sys.argv[2])
        test_rounds = int(sys.argv[3])
        test_accuracy(test_number, test_rounds)
    elif sys.argv[1] == '-s':
        test_number = int(sys.argv[2])
        test_statistics(test_number)
    elif sys.argv[1] == '-speed':
        test_number = int(sys.argv[2])
        test_speed(test_number)

# 179425457
# 141642509
# 92171017
# 561, 8911 (числа Кармайкла)
# 3277 - полупростое число = p * q (простые)
# 349 181 521 761 907 911
