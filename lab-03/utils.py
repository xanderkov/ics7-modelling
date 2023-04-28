from math import *
from numpy import arange
from table_reader import TableReader
from newton_poly import newton_interpol

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
table_l = TableReader("f1.txt")
newton_table1 = table_l.make_table_from_file()
list_lam = table_l.dot_arr

table_k = TableReader("f2.txt")
newton_table2 = table_k.make_table_from_file()
list_k = table_k.dot_arr

# print_newt_table(newton_table1)


def integration(ys, z):
    s = 0.5 * (k_f(ys[0]) * z[0] * (ys[0] ** 4 - T0 ** 4) + k_f(ys[N]) * z[N] * (ys[N] ** 4 - T0 ** 4))
    for n, z_n in enumerate(z[1:N]):
        s += z_n * k_f(ys[n]) * (ys[n] ** 4 - T0 ** 4)
    return h * s


# ============================ = Краевые условия = =================================


def z(n):
    return r0 / R + h * n


def z_m12(n, z_0):
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
        B[n] = (z_n_p_half * kappa_n_p_half + z_n_m_half * kappa_n_m_half) / (
                    R ** 2 * h) + 2 * np ** 2 * sigma * k_t(
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


T_sp = [T_start] * (N + 1)
z_sp = arange(r0 / R, 1 + EPS, h)

z_0 = z(0)
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
_t = z_m12(n, z_0)
_t1 = kappa_m12(n)

