from typing import Callable
from models import Node, Numerical
import numpy as np


def solve_euler(node: Node, x: float, y: float, f: Callable[[float, float], float]):
    answer = []
    n = len(np.arange(node.x_min, node.x_max, node.h))
    for i in range(n):
        try:
            y += node.h * f(x, y)
            answer.append(y)
            x += node.h
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