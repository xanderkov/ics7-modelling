from typing import Callable
from dataclasses import dataclass


@dataclass
class Node:
    x_min: float
    x_max: float
    h: float


@dataclass
class Numerical:
    f_eyler: Callable[[float, float], float]
    approx_picar_1: Callable[[float, float], float]
    approx_picar_2: Callable[[float], float]
    approx_picar_3: Callable[[float], float]
    approx_picar_4: Callable[[float], float]