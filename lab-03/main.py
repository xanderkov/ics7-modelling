from utils import *
from numpy import arange
import matplotlib.pyplot as plt

step = 1e-6


def img_routine(plot, name):
    plot.legend()
    plot.grid()
    plt.ylabel(name)
    plt.xlabel("z, б/р")

    plt.show()


def graph(num=1, spt = None, spi = None):
    num = str(num)
    s1 = 300
    s2 = 2400
    n = 100
    h = (s2 - s1 )/ n
    spt = arange(300, 2400 + h, h)
    spi = []
    for i in spt:
        t = newton_interpol(newton_table2, log(i), 2)
        spi.append(exp(t))
    fig = plt.figure(figsize=(10, 7))
    plot = fig.add_subplot()
    plot.plot(spt, spi,  label = "runge 2", c = 'r')
    plot.plot([i[2] for i in list_k], [i[3] for i in list_k], '*')
    img_routine(plot, "Сила тока I, А")
    return


def graph1():
    fig = plt.figure(figsize=(10, 7))
    plot = fig.add_subplot()
    plot.plot(z_sp, T_sp,  label = "T", c = 'r')
    img_routine(plot, "Температура T, K")
    return


def main():
    graph1()


if __name__ == "__main__":
    main()
