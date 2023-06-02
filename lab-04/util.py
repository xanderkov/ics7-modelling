from math import *
from numpy import arange
import copy


V = 2
EPS = 1e-6
EPS2 = 1e-4
EPS1 = 0.05
T_start = 300
N = 50


np = 1.4
r0 = 0.35
R = 0.5
T0 = 300
sigma = 5.668e-12
F0 = 100
Fmax = 5000
tmax = 150
alpha = 0.05

a2 = 2.049
b2 = 0.563e-3
c2 = 0.528e5
m2 = 1

h = (1-r0/R) / N
cons = np*np*sigma*h


def print_newt_table(table):
    n = len(table)
    for i in range(n):
        for j in range(n - i + 1):
            print("{:<15.8f}".format(table[i][j]), end = ' ')
        print()

def load_dots_from_file(filename, ind = 1):
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
                dots.append([x, y])
            except:
                print("File wasn't properly read")
                break
            line = f.readline()
    f.close()
    return dots

def make_newton_table(dot_arr):
    n = len(dot_arr)
    m = len(dot_arr) + 1
    table  = [0] * n
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

        for i in dot_list[last - before + 1 : last + after + 1]:
            arr.append(i)
        return arr


def newton_polinome(table, x):
    y = 0
    multiplier = 1
    for i in range(1, len(table[0])):
        y += multiplier * table[0][i]
        multiplier *= x - table[i - 1][0]
    return y


def newton_interpol(input_list, x, n):
    arr = choose_dots(x, input_list, n)
    if len(arr) != n:
        return None, None
    newton_table = make_newton_table(arr)

    return newton_polinome(newton_table, x)

def integration(ys, z):
    s = 0.5 * (k_f(ys[0]) * z[0] * (ys[0] ** 4 - T0 ** 4) + k_f(ys[N]) * z[N] * (ys[N] ** 4 - T0 ** 4))
    for n, z_n in enumerate(z[1:N]):
        s += z_n * k_f(ys[n]) * (ys[n] ** 4 - T0 ** 4)
    return h * s

# ============================ = Функции = =================================
def z(n):
    return r0/R + h*n
z_0 = z(0)


def lambda_t(y):
    t = newton_interpol(list_lam, y, 2)
    return t
def k_t(y):
    t = newton_interpol(list_k, y, 2)
    return t

def k_f(y):
    t = newton_interpol(list_k, y, 2)
    return t

def c(T): 
    return a2 + b2*T**m2 - c2/T/T

def F(t):
    return Fmax/tmax * t * exp(1 - t/tmax)
    # return 300




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
    A[N] -= np ** 2 * sigma * h * (k_t(ys[N - 1]) + k_t(ys[N])) / 2 * ((ys[N - 1] + ys[N]) / 2) ** 3 * (z[N - 1] + z[N]) / 2 / 2
    return A

def get_B(z, ys):
    B = [0] * (N + 1)
    h = z[1] - z[0]
    z_n_p_half = (z[1] + z[0]) / 2
    kappa_n_p_half = (lambda_t(ys[1]) + lambda_t(ys[0])) / 2
    B[0] = z_n_p_half * kappa_n_p_half / (R ** 2 * h) + np ** 2 * sigma * h * ((k_t(ys[0]) + k_t(ys[1])) / 2 * z_n_p_half * ((ys[0] + ys[1]) / 2) ** 3 / 2 + k_t(ys[0]) * z[0] * ys[0] ** 3)
    for n in range(1, N):
        z_n_m_half = (z[n - 1] + z[n]) / 2
        kappa_n_m_half = (lambda_t(ys[n - 1]) + lambda_t(ys[n])) / 2
        z_n_p_half = (z[n + 1] + z[n]) / 2
        kappa_n_p_half = (lambda_t(ys[n + 1]) + lambda_t(ys[n])) / 2
        B[n] = (z_n_p_half * kappa_n_p_half + z_n_m_half * kappa_n_m_half) / (R ** 2 * h) + 2 * np ** 2 * sigma * k_t(ys[n]) * ys[n] ** 3 * ((z_n_p_half) ** 2 - (z_n_m_half) ** 2)
    z_n_m_half = (z[N - 1] + z[N]) / 2
    kappa_n_m_half = (lambda_t(ys[N - 1]) + lambda_t(ys[N])) / 2
    B[N] = z[N] * alpha / R + z_n_m_half * kappa_n_m_half / (R ** 2 * h) + np ** 2 * sigma * h * ((k_t(ys[N]) + k_t(ys[N - 1])) / 2 * z_n_m_half * ((ys[N] + ys[N - 1]) / 2) ** 3 / 2 + k_t(ys[N]) * z[N] * ys[N] ** 3)
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

def get_D(z, ys, t):
    D = [0] * (N + 1)
    h = z[1] - z[0]
    D[0] = z[0] * F(t) / R + np ** 2 * sigma * h * T0 ** 4 * (k_t(ys[0]) * z[0] + (k_t(ys[0]) + k_t(ys[1])) / 2 * (z[0] + z[1]) / 2)
    for n in range(1, N):
        z_n_m_half = (z[n - 1] + z[n]) / 2
        z_n_p_half = (z[n + 1] + z[n]) / 2
        D[n] = 2 * np ** 2 * sigma * T0 ** 4 * k_t(ys[n]) * ((z_n_p_half) ** 2 - (z_n_m_half) ** 2)
    D[N] = np ** 2 * sigma * h * T0 ** 4 * (k_t(ys[N]) * z[N] + (k_t(ys[N]) + k_t(ys[N - 1])) / 2 * (z[N] + z[N - 1]) / 2) + z[N] * alpha * T0 / R
    return D

def get_A2(z, ys, ynext, tau):
    A = get_A(z, ynext)
    for n in range(1, N):
        A[n] *= tau
    A[N] = A[N]*tau + h/4*(c(ynext[N]) + c(ynext[N-1]))/2
    # Раскоментить, если T в левом крае константа
    A[N] = 0
    return A

def get_B2(z, ys, ynext, tau):
    B = get_B(z, ynext)
    for n in range(1, N):
        B[n] = B[n] * tau + c(ynext[n]) * h
    B[0] = -B[0] * tau - h/2*(c(ynext[0]) + (c(ynext[1]) + c(ynext[0]))/4)
    B[0] *= -1
    B[N] = -B[N]*tau - h/2*(c(ynext[N]) + (c(ynext[N]) + c(ynext[N-1]))/4)
    B[N] *= -1
    #раскомментить, если нужно задание 5
    # B[0] = -1
    B[N] = -1
    return B

def get_E2(z, ys, ynext, tau):
    E = get_C(z, ynext)
    for n in range(1, N):
        E[n] *= tau
    E[0] = E[0] * tau - h/4 * ((c(ynext[1]) + c(ynext[0]))/2)
    #раскомментить, если нужно задание 5
    # E[0] = 0
    return E

def get_D2(z, ys,ynext,tau, t):
    D = get_D(z, ynext, t)
    for n in range(1, N):
        D[n] = D[n] * tau + c(ynext[n]) * ys[n] * h
    D[0] = -D[0] * tau - h/2*(c(ynext[0])*ys[0] + (c(ynext[1]) + c(ynext[0]))/2 * (ys[0] + ys[1])/2)
    D[0] *=-1
    D[N] = -D[N] * tau - h/2*(c(ynext[N])*ys[N] + (c(ynext[N]) + c(ynext[N-1]))/2 * (ys[N] + ys[N-1])/2)
    D[N] *=-1
    #раскомментить, если нужно задание 5
    # D[0] = -2000
    D[N] = -3300
    return D

def progon(z, ys, ynext, tau, t):
    A = get_A2(z, ys, ynext, tau)
    B = get_B2(z, ys, ynext, tau)
    C = get_E2(z, ys, ynext, tau)
    D = get_D2(z, ys, ynext, tau, t)
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
    return n == N + 1
