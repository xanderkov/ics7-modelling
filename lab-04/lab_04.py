from random import randint
from math import *
from util import *
from numpy import arange
import pandas as pd
import prettytable
import matplotlib.pyplot as plt

# Раскомментить, если нужно задание 5
# M = 3000 
M = 100

def img_routine(plot, name, name1):
    plot.legend()
    plot.grid()
    plt.ylabel(name)
    plt.xlabel(name1)

    plt.show()

def graph2(z_sp, res, tm, tau, m):
    fig = plt.figure(figsize=(10, 7)) 
    plot = fig.add_subplot()
    sp = [tau*i for i in range(int(m))]
    sp2 = [i[0] for i in res]
    plot.plot(sp, sp2,  label = "T[0]", c = "m")
    img_routine(plot, "Температура T[0], K","t, мкс")
    return

def graph1(z_sp, res, tau):
    n = 10
    h = (M) // n
    fig = plt.figure(figsize=(10, 7)) 
    plot = fig.add_subplot()
    for k in range(n):
        plot.plot(z_sp, res[1 + h * k],  label = "t = " + str((1 + h * k) * tau))
    img_routine(plot, "Температура T, K", "z, б/р")
    return

def write_res(z_sp, res, tm, tau,m):
    a = input("Введите название файла для сохранения данных:")
    f = open(a, 'w')
    f.write(str(tm) + ' ' + str(tau) + ' '+str(m)+'\n')
    st = ''
    for i in z_sp:
        st += str(i) + ' '
    f.write(st + '\n')
    for line in res:
        st = ''
        for i in line:
            st += str(i) + ' '
        f.write(st + '\n')
    f.close()

def read_res(f):
    tm, tau,m = list(map(float, f.readline().split()))
    z_sp = list(map(float, f.readline().split()))
    res = []
    for line in f.readlines():
        s = list(map(float, line.split()))
        res.append(s)
    return z_sp, res,tm,tau,m

def main():
    # Раскоментить, если нужно задание 5
    # tlimit = tmax*20 
    tlimit = tmax * 1.5
    tau = tlimit/M 
    t = 0
    m = 0
    T_sp = [T_start] * (N + 1)
    res = []
    z_sp = arange(r0/R, 1+EPS, h)
    while m < M:
        res.append(T_sp)
        print("m = ",m)
        m+=1
        t+=tau
        run = True
        T_sp_next = copy.deepcopy(T_sp)
        cnt = 1
        while run:
            cnt+=1
            tmp = progon(z_sp,T_sp, T_sp_next, tau, t%1000)
            run = not check(T_sp_next, tmp)
            T_sp_next = [(tmp[i] + T_sp_next[i]) / 2 for i in range(N+1)]
        
        T_sp = T_sp_next
    write_res(z_sp, res, tlimit, tau, M)
    return z_sp, res,tlimit, tau,M
    

    
#В файле "1" графики для 3,4 задания;
#В файле "5" графики для 5 задания 
#Ввод названия несуществующего файла запускает вычисления
a = input("Введите название файла:")
f = None
try:
    f = open(a, 'r')
    z_sp, res,tm, tau,m = read_res(f)
    f.close()
except FileNotFoundError:
    z_sp, res,tm, tau,m = main()

graph1(z_sp, res, tau)
graph2(z_sp, res, tm, tau, m)
pass