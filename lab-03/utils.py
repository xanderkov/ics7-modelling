from math import *
from numpy import arange


EPS = 1e-6
EPS1 = 0.05
T_start = 300
N = 1000
# N = 32 - minimum, при котором сходится

np = 1.4
r0 = 0.35
R = 0.5
T0 = 300
sigma = 5.668e-12
F0 = 100
# F0 = -7 min, отрицательный -- положительная производная
alpha = 0.05
# при умножении на 3 температура всего цилиндра снижается

h = (1 - r0 / R) / N
cons = np * np * sigma * h


# h = 1e-6
# N = int((R-r0)/h)


def print_newt_table(table):
    n = len(table)
    for i in range(n):
        for j in range(n - i + 1):
            print("{:<15.8f}".format(table[i][j]), end=' ')
        print()


def load_dots_from_file(filename, ind=1):
    try:
        f = open(filename)
    except:
        print("File doesn't exist")
        return []
    dots = []
    line = f.readline()
    if ind == 1:
        while line:
            try:
                x, y = map(float, line.split())
                dots.append([log(x), log(y), x, y])
                # dots.append([x, y, x, y])
            except:
                print("File wasn't properly read")
                break
            line = f.readline()
    f.close()
    return dots


def make_newton_table(dot_arr):
    n = len(dot_arr)
    m = len(dot_arr) + 1
    table = [0] * n
    for i in range(n):
        table[i] = [0] * m
    for i in range(n):
        table[i][0] = dot_arr[i][0]
        table[i][1] = dot_arr[i][1]
    for j in range(2, n + 1):
        for i in range(n - j + 1):
            table[i][j] = (table[i][j - 1] - table[i + 1][j - 1]) / (table[i][0] - table[i + (j - 1)][0])
    return table


list_lam = load_dots_from_file("f1.txt")
list_k = load_dots_from_file("f2.txt")

newton_table1 = make_newton_table(list_lam)
newton_table2 = make_newton_table(list_k)


# print_newt_table(newton_table1)


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
    arr = choose_dots(x, input_list, n)
    if len(arr) != n:
        return None, None
    newton_table = make_newton_table(arr)

    # print_newt_table(newton_table)
    return newton_polinome(newton_table, x)


def integration(ys, z):
    s = 0.5 * (k_f(ys[0]) * z[0] * (ys[0] ** 4 - T0 ** 4) + k_f(ys[N]) * z[N] * (ys[N] ** 4 - T0 ** 4))
    for n, z_n in enumerate(z[1:N]):
        s += z_n * k_f(ys[n]) * (ys[n] ** 4 - T0 ** 4)
    return h * s


# ============================ = Краевые условия = =================================
def z(n):
    return r0 / R + h * n


z_0 = z(0)


def z_p12(n):
    return z_0 + (n + 0.5) * h


def z_m12(n):
    return z_0 + (n - 0.5) * h


def kappa_p12(n):
    return lamda_p12(n)


def kappa_m12(n):
    return lamda_m12(n)


def lamda_p12(n):
    return (lamda(n) + lamda(n + 1)) / 2


def lamda_m12(n):
    return (lamda(n) + lamda(n - 1)) / 2


def k_p12(n):
    return (k(n) + k(n + 1)) / 2


def k_m12(n):
    return (k(n) + k(n - 1)) / 2


def y_p12(n):
    return (y(n) + y(n + 1)) / 2


def y_m12(n):
    return (y(n) + y(n - 1)) / 2


def y(n):
    return T_sp[n]


def lamda(n):
    T = T_sp[n]
    t = newton_interpol(list_lam, log(T), 3)
    return exp(t)


def k(n):
    T = T_sp[n]
    t = newton_interpol(list_k, log(T), 3)
    return exp(t)


def lambda_t(y):
    t = newton_interpol(list_lam, log(y), 3)
    return exp(t)


def k_t(y):
    t = newton_interpol(list_k, log(y), 3)
    return exp(t)


def k_f(y):
    t = newton_interpol(list_k, log(y), 3)
    return exp(t)


T_sp = [T_start] * (N + 1)
z_sp = arange(r0 / R, 1 + EPS, h)

z_12 = z_0 + h / 2
kappa_12 = kappa_p12(0)
k_0 = k(0)
k_12 = k_p12(0)
y_12 = y_p12(0)
M0 = -z_12 * kappa_12 / R / R / h - cons * (k_12 * z_12 * y_12 ** 3 / 2 + k_0 * z_0 * y(0))
K0 = z_12 * kappa_12 / R / R / h - cons * k_12 * z_12 * y_12 / 2
P0 = -(z_0 * F0 / R + cons * T0 ** 4 * (k_0 * z_0 + k_12 * z_12))

z_N = z(N)
z_N_m12 = z_N - h / 2
kappa_N_m12 = kappa_m12(N)
k_N = k(N)
k_N_m12 = k_m12(N)
y_N_m12 = y_m12(N)
y_N = y(N)
M_N = 1 / R / R / h * z_N_m12 * kappa_N_m12 - cons * k_N_m12 * y_N_m12 ** 3 * z_N_m12 / 2
K_N = -1 / R * z_N * alpha - 1 / R / R / h * z_N_m12 * kappa_N_m12 - cons * (
            k_N_m12 * y_N_m12 ** 3 / 2 + k_N * y_N ** 3 * z_N)
P_N = -cons * T0 ** 4 * (k_N * z_N + k_N_m12 * z_N_m12) - z_N * alpha * T0 / R


# ============================== = Метод прогонки = ===============================================
# A = [0] * (N + 1)
# B = [0] * (N + 1)
# C = [0] * (N + 1)
# D = [0] * (N + 1)

def get_A(z, ys):
    A = [0] * (N + 1)
    h = z[1] - z[0]
    for n in range(1, N + 1):
        z_n_m_half = (z[n - 1] + z[n]) / 2
        kappa_n_m_half = (lambda_t(ys[n - 1]) + lambda_t(ys[n])) / 2
        A[n] = z_n_m_half * kappa_n_m_half / (R ** 2 * h)
    A[N] -= np ** 2 * sigma * h * (k_t(ys[N - 1]) + k_t(ys[N])) / 2 * ((ys[N - 1] + ys[N]) / 2) ** 3 * (
                z[N - 1] + z[N]) / 2 / 2
    return A


def get_B(z, ys):
    B = [0] * (N + 1)
    h = z[1] - z[0]
    z_n_p_half = (z[1] + z[0]) / 2
    kappa_n_p_half = (lambda_t(ys[1]) + lambda_t(ys[0])) / 2
    B[0] = z_n_p_half * kappa_n_p_half / (R ** 2 * h) + np ** 2 * sigma * h * (
                (k_t(ys[0]) + k_t(ys[1])) / 2 * z_n_p_half * ((ys[0] + ys[1]) / 2) ** 3 / 2 + k_t(ys[0]) * z[0] * ys[
            0] ** 3)
    for n in range(1, N):
        z_n_m_half = (z[n - 1] + z[n]) / 2
        kappa_n_m_half = (lambda_t(ys[n - 1]) + lambda_t(ys[n])) / 2
        z_n_p_half = (z[n + 1] + z[n]) / 2
        kappa_n_p_half = (lambda_t(ys[n + 1]) + lambda_t(ys[n])) / 2
        B[n] = (z_n_p_half * kappa_n_p_half + z_n_m_half * kappa_n_m_half) / (R ** 2 * h) + 2 * np ** 2 * sigma * k_t(
            ys[n]) * ys[n] ** 3 * ((z_n_p_half) ** 2 - (z_n_m_half) ** 2)
    z_n_m_half = (z[N - 1] + z[N]) / 2
    kappa_n_m_half = (lambda_t(ys[N - 1]) + lambda_t(ys[N])) / 2
    B[N] = z[N] * alpha / R + z_n_m_half * kappa_n_m_half / (R ** 2 * h) + np ** 2 * sigma * h * (
                (k_t(ys[N]) + k_t(ys[N - 1])) / 2 * z_n_m_half * ((ys[N] + ys[N - 1]) / 2) ** 3 / 2 + k_t(ys[N]) * z[
            N] * ys[N] ** 3)
    return B


def get_C(z, ys):
    C = [0] * (N + 1)
    h = z[1] - z[0]
    for n in range(N):
        z_n_p_half = (z[n + 1] + z[n]) / 2
        kappa_n_p_half = (lambda_t(ys[n + 1]) + lambda_t(ys[n])) / 2
        C[n] = z_n_p_half * kappa_n_p_half / (R ** 2 * h)
    C[0] -= np ** 2 * sigma * h * (k_t(ys[1]) + k_t(ys[0])) / 2 * ((ys[1] + ys[0]) / 2) ** 3 * (z[1] + z[0]) / 2 / 2
    return C


def get_D(z, ys):
    D = [0] * (N + 1)
    h = z[1] - z[0]
    D[0] = z[0] * F0 / R + np ** 2 * sigma * h * T0 ** 4 * (
                k_t(ys[0]) * z[0] + (k_t(ys[0]) + k_t(ys[1])) / 2 * (z[0] + z[1]) / 2)
    for n in range(1, N):
        z_n_m_half = (z[n - 1] + z[n]) / 2
        z_n_p_half = (z[n + 1] + z[n]) / 2
        D[n] = 2 * np ** 2 * sigma * T0 ** 4 * k_t(ys[n]) * ((z_n_p_half) ** 2 - (z_n_m_half) ** 2)
    D[N] = np ** 2 * sigma * h * T0 ** 4 * (
                k_t(ys[N]) * z[N] + (k_t(ys[N]) + k_t(ys[N - 1])) / 2 * (z[N] + z[N - 1]) / 2) + z[N] * alpha * T0 / R
    return D


def progon(z, ys):
    A = get_A(z, ys)
    B = get_B(z, ys)
    C = get_C(z, ys)
    D = get_D(z, ys)
    ksi = [0] * (N + 1)
    eta = [0] * (N + 1)
    ys1 = [T_start] * (N + 1)
    for n in range(N):
        low = B[n] - A[n] * ksi[n]
        ksi[n + 1] = C[n] / low
        eta[n + 1] = (D[n] + A[n] * eta[n]) / low

    ys1[N] = (A[N] * eta[N] + D[N]) / (B[N] - A[N] * ksi[N])
    for n in range(N, 0, -1):
        ys1[n - 1] = ksi[n] * ys1[n] + eta[n]
    return ys1


def check(ys, ys1):
    n = 0
    while n < N + 1 and abs((ys[n] - ys1[n]) / ys[n]) < EPS:
        n += 1
    if n == N + 1:
        f1 = r0 * F0 - R * alpha * (ys[N] - T0)
        f2 = 4 * np * np * sigma * R ** 2 * integration(ys, z_sp)
        print(f"f1 = {f1}, f2 = {f2}")
        if abs(f1) < EPS or abs((f1 - f2) / f1) < EPS1:
            return True
    return False


run = True
T_sp = progon(z_sp, T_sp)
cnt = 1
while run:
    cnt += 1
    tmp = progon(z_sp, T_sp)
    run = not check(T_sp, tmp)
    T_sp = tmp
    print("cnt = ", cnt)
n = 1
_t = z_m12(n)
_t1 = kappa_m12(n)