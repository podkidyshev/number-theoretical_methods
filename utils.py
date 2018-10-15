import euclid

MAX_B_SIZE = 5000
MAX_SUBSET = 30


def get_inverse(a, m):
    if euclid.euclid(a, m) != 1:
        raise ValueError('Не существует обратного элемента для a={} по модулю m={}'.format(a, m))
    d, x, y = euclid.eeuclid(a, m)
    assert d == 1
    return x % m


def fac2k(a):
    k = 0
    while a & 1 == 0:
        a >>= 1
        k += 1
    return a, k


def legendre(a, n):
    a %= n
    if a == 0:
        return 0
    elif a == 1:
        return 1
    return pow(a, (n - 1) // 2, n)


def jacobi(a, n, g=1):
    """Маховенко Е.Б. Теоретико-числовые методы в криптографии, стр 61-62"""
    if a == 0:
        return 0
    elif a == 1:
        return g

    a1, k = fac2k(a)

    if k & 1 == 0 or n % 8 == 1 or n % 8 == 7:
        s = 1
    else:
        s = -1

    if a1 == 1:
        return g * s

    if n % 4 == 3 and a1 % 4 == 3:
        s = -s

    return jacobi(n % a1, a1, g * s)


def generate_base(desired_count):
    base = [2, 3, 5]
    counter = 7
    while len(base) < desired_count:
        for b in base:
            if counter % b == 0:
                break
            if b > counter // 2:
                base.append(counter)
                break
        counter += 1
    return base


def matrix_print(a):
    for row in a:
        print(row)


if __name__ == '__main__':
    pass
