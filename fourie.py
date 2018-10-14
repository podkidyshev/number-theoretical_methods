def find_root(p, deg):
    for w in range(2, p):
        if pow(w, deg, p) == 1:
            return w
    else:
        raise ArithmeticError('Нет корня степени {} из {}'.format(deg, p))


if __name__ == '__main__':
    p_test = 1300633
    deg_test = p_test - 1
    a = find_root(p_test, deg_test)
    assert pow(a, deg_test, p_test) == 1
    print(a)
