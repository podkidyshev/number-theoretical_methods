import copy
from itertools import chain, combinations


def zero_row(row):
    return all(el == 0 for el in row)


def sum_rows(*rows):
    # результат будет в первом аргументе (для удобства он будет возвращен)
    for row in rows[1:]:
        for idx in range(len(row)):
            rows[0][idx] = (rows[0][idx] + row[idx]) % 2
    return rows[0]


def all_subsets(elements, *args):
    rng = range(*args) if len(args) else range(0, len(elements) + 1)
    return chain(*map(lambda x: combinations(elements, x), rng))


def gen_gaussian_rec(es, col_num, row_num, stack, zeros):
    # 1. поксорить
    indices = set(range(len(es))) - stack

    row_set = []
    for row_idx in indices:
        if es[row_idx][col_num] == 1 and all(el == 0 for el in es[row_idx][:col_num]):
            row_set.append(row_idx)

    for row_subset in all_subsets(row_set, 1, len(row_set) + 1, 2):
        sum_rows(es[row_num], *[es[r] for r in row_subset])
        stack = stack.union(row_subset)

        # 2. выдать ответы если есть
        if all(el == 0 for el in es[row_num]):
            for zero_subset in all_subsets(zeros):
                yield tuple(stack) + zero_subset

        if any(el == 1 for el in es[row_num]):
            col_num_rec = col_num
            while es[row_num][col_num_rec] == 0:
                col_num_rec += 1

            yield from gen_gaussian_rec(es, col_num_rec, row_num, stack, zeros)

        sum_rows(es[row_num], *[es[r] for r in row_subset])
        stack = stack.difference(row_subset)


def gen_gaussian(es_r):
    es = copy.deepcopy(es_r)
    numbers = set()
    zeros = list(filter(lambda rr: zero_row(es[rr]), range(len(es))))
    # композиции нулевых строк
    for subset in all_subsets(zeros, 1, len(zeros) + 1):
        yield subset
    # комбинации нулевых и ненулевых решений
    for idx, row in enumerate(es):
        if idx in zeros:
            continue
        idx_1 = 0
        while row[idx_1] == 0:
            idx_1 += 1
        numbers.add(idx)
        yield from gen_gaussian_rec(es, idx_1, idx, numbers, zeros)
        numbers.remove(idx)
