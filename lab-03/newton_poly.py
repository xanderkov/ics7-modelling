from table_reader import TableReader
# ====================================================== == Ньютоновская интерполяция ================================================= ===================


def choose_dots(x, dot_list, n):
    if n > len(dot_list):
        return []
    else:
        before = n // 2
        after = n - before
        arr = []
        last = 0
        for i in range(len(dot_list)):
            if dot_list[i][0] < x:
                last = i
            else:
                break

        if (last + 1) < before:
            before = last + 1
            after = n - before
        elif (len(dot_list) - last - 1) < after:
            after = len(dot_list) - last - 1
            before = n - after

        for i in dot_list[last - before + 1: last + after + 1]:
            arr.append(i)
        return arr


def newton_polinome(table, x):
    y = 0
    multiplier = 1
    for i in range(1, len(table[0])):
        y += multiplier * table[0][i]
        # print(y)
        multiplier *= x - table[i - 1][0]
    return y


def newton_interpol(input_list, x, n):
    # sorted_list = sorted(input_list, key = lambda x: x[0])
    dot_arr = choose_dots(x, input_list, n)
    if len(dot_arr) != n:
        return None, None
    table = TableReader(arr=dot_arr)
    newton_table = table.make_newton_table()

    # print_newt_table(newton_table)
    return newton_polinome(newton_table, x)
