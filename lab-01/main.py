from prettytable import PrettyTable
import numpy as np
import math


def euler(n, h, x, y, f):
    answer = []

    for i in range(n):
        try:
            y += h * f(x, y)
            answer.append(y)
            x += h
        except OverflowError:
            answer.append("Over")

    return answer


def f_1(x, u):
	return x + u**2


def first_get_anal_solve(x_min, x_max, h) -> list:
    fx = lambda x: 3 * np.exp(x) - x**2 - 2*x -2
    anal = []
    for i in np.arange(x_min, x_max, h):
        anal.append(fx(i))
    return anal

def get_eyler_solve(x_min, x_max, h) -> list:
    eyler = []
    fx = lambda x: x**2
    for i in np.arange(x_min, x_max, h):
        eyler.append(fx(i))
    return eyler

def get_pikar_solve_1(x_min, x_max, h) -> list:
    pikar_1 = []
    fx = lambda x: x**2
    for i in np.arange(x_min, x_max, h):
        pikar_1.append(fx(i))
    return pikar_1

def get_pikar_solve_2(x_min, x_max, h) -> list:
    pikar_2 = []
    fx = lambda x: x**2
    for i in np.arange(x_min, x_max, h):
        pikar_2.append(fx(i))
    return pikar_2

def get_pikar_solve_3(x_min, x_max, h) -> list:
    pikar_3 = []
    fx = lambda x: x**2
    for i in np.arange(x_min, x_max, h):
        pikar_3.append(fx(i))
    return pikar_3

def get_pikar_solve_4(x_min, x_max, h) -> list:
    pikar_4 = []
    fx = lambda x: x**2
    for i in np.arange(x_min, x_max, h):
        pikar_4.append(fx(i))
    return pikar_4


def solve_first_task(x_min: float, x_max: float, h: float) -> PrettyTable:
    task = PrettyTable()
    
    task.field_names = ["Аргумент", "Аналит.", "Эйлер", "Пикар-1", "Пикар-2", "Пикар-3", "Пикар-4"]
    
    anal = first_get_anal_solve(x_min, x_max, h)
    eyler = get_eyler_solve(x_min, x_max, h)
    pikar_1 = get_pikar_solve_1(x_min, x_max, h)
    pikar_2 = get_pikar_solve_2(x_min, x_max, h)
    pikar_3 = get_pikar_solve_3(x_min, x_max, h)
    pikar_4 = get_pikar_solve_4(x_min, x_max, h)  
    
    for i in range(len(anal)):
        task.add_row([i, anal[i], eyler[i], pikar_1[i], pikar_2[i], pikar_3[i], pikar_4[i]])
    
    return task


def main():
    x_min, x_max, h = -1, 1, 0.1
    # x_min, x_max, h = input("Max, Min, Step: ").split()  
    
    first_task = solve_first_task(x_min, x_max, h)
    print(first_task)



if __name__ == "__main__":
    main()
