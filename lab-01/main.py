from prettytable import PrettyTable
import numpy as np
import math
from typing import Callable
from dataclasses import dataclass


@dataclass
class Node:
    x_min: float
    x_max: float
    h: float



def solve_euler(node: Node, x0: float, y0: float, f: Callable[[float, float], float]) -> list:
    answer = []
    for x0 in np.arange(node.x_min, node.x_max, node.h):
        try:
            y0 += node.h * f(x0, y0)
            answer.append(y0)
        except:
            answer.append("Переполнение")

    return answer


def get_anal_solve(node: Node, f: Callable[[float], float]) -> list:
    anal = []
    for i in np.arange(node.x_min, node.x_max, node.h):
        anal.append(f(i))
    return anal


def solve_picar(node: Node, aprox_func: Callable[[float], float]) -> list:
    picar = []
    u = 0
    for x in np.arange(0, np.abs(node.x_max), np.abs(node.h)):
        picar.append(u)
        u = aprox_func(x)
    return picar


def solve_first_task(node: Node) -> PrettyTable:
    task = PrettyTable()
    
    f_1 = lambda x, u: x + np.power(u, 2)
    f_anal = lambda x: 3 * np.exp(x) - x**2 - 2*x -2
    approx_picar_1 = lambda u: 1 + u + np.power(u, 3) / 3
    
    task.field_names = ["Аргумент", "Аналит.", "Эйлер", "Пикар-1", "Пикар-2", "Пикар-3", "Пикар-4"]
    
    anal = get_anal_solve(node, f_anal)
    eyler = solve_euler(node, 0, 1, f_1)
    pikar_1 = solve_picar(node, approx_picar_1)
    pikar_2 = solve_picar(node, approx_picar_1)
    pikar_3 = solve_picar(node, approx_picar_1)
    pikar_4 = solve_picar(node, approx_picar_1)
    
    for i in range(len(anal)):
        row = [i, anal[i], eyler[i]]
        try:
            picar = [pikar_1[i], pikar_2[i], pikar_3[i], pikar_4[i]]
        except:
            picar = ["Ошибка", "Ошибка","Ошибка","Ошибка"]
        row.extend(picar)
        task.add_row(row)
    
    return task


def main():
    x_min, x_max, h = -1, 1, 0.1
    # x_min, x_max, h = input("Max, Min, Step: ").split()  
    
    node = Node(x_min=x_min, x_max=x_max, h=h)
    
    first_task = solve_first_task(node)
    print(first_task)



if __name__ == "__main__":
    main()
