from prettytable import PrettyTable
import numpy as np
from models import Node, Numerical
from math_methods import solve_euler, solve_picar, get_anal_solve


def get_numerical_methods(node: Node, Numerical: Numerical, x0: float, y0: float):
    eyler = solve_euler(node, x0, y0, Numerical.f_eyler)
    pikar_1 = solve_picar(node, Numerical.approx_picar_1)
    pikar_2 = solve_picar(node, Numerical.approx_picar_2)
    pikar_3 = solve_picar(node, Numerical.approx_picar_3)
    pikar_4 = solve_picar(node, Numerical.approx_picar_4)
    x = np.arange(node.x_min, node.x_max, node.h)

    return eyler, pikar_1, pikar_2, pikar_3, pikar_4, x


def solve_first_task(node: Node) -> PrettyTable:
    task = PrettyTable()

    f_eyler = lambda x, u: (x + u**2)
    f_anal = lambda x: 3 * np.exp(x) - x**2 - 2 * x - 2
    approx_picar_1 = lambda u: 1 + u + np.power(u, 3) / 3
    approx_picar_2 = (
        lambda u: approx_picar_1(u) + np.power(u, 2) / 2 + np.power(u, 4) / 12
    )
    approx_picar_3 = (
        lambda u: approx_picar_2(u) + np.power(u, 3) / 6 + np.power(u, 5) / 60
    )
    approx_picar_4 = (
        lambda u: approx_picar_3(u) + np.power(u, 4) / 25 + np.power(u, 6) / 360
    )

    numerical = Numerical(
        f_eyler=f_eyler,
        approx_picar_1=approx_picar_1,
        approx_picar_2=approx_picar_2,
        approx_picar_3=approx_picar_3,
        approx_picar_4=approx_picar_4,
    )

    task.field_names = [
        "Аргумент",
        "Аналит.",
        "Эйлер",
        "Пикар-1",
        "Пикар-2",
        "Пикар-3",
        "Пикар-4",
    ]

    anal = get_anal_solve(node, f_anal)
    eyler, pikar_1, pikar_2, pikar_3, pikar_4, x = get_numerical_methods(
        node, numerical, 1, 0
    )

    for i in range(len(anal)):
        row = [x[i], anal[i], eyler[i]]
        try:
            picar = [pikar_1[i], pikar_2[i], pikar_3[i], pikar_4[i]]
        except:
            picar = ["Ошибка", "Ошибка", "Ошибка", "Ошибка"]
        row.extend(picar)
        task.add_row(row)

    return task


def solve_second_task(node: Node) -> PrettyTable:
    task = PrettyTable()

    f_eyler = lambda x, u: (u**3 + 2*x*u)
    f_anal = lambda x: np.exp(x**2) - (x**2 + 1) / 2
    approx_picar_1 = lambda u: 0.5 + np.power(u, 2) / 2 + np.power(u, 4) / 4
    approx_picar_2 = (
        lambda u: approx_picar_1(u) + np.power(u, 4) / 4 + np.power(u, 6) / 12
    )
    approx_picar_3 = (
        lambda u: approx_picar_2(u) + np.power(u, 6) / 12 + np.power(u, 8) / 48
    )
    approx_picar_4 = (
        lambda u: approx_picar_3(u) + np.power(u, 8) / 48 + np.power(u, 10) / 240
    )

    numerical = Numerical(
        f_eyler=f_eyler,
        approx_picar_1=approx_picar_1,
        approx_picar_2=approx_picar_2,
        approx_picar_3=approx_picar_3,
        approx_picar_4=approx_picar_4,
    )

    task.field_names = [
        "Аргумент",
        "Аналит.",
        "Эйлер",
        "Пикар-1",
        "Пикар-2",
        "Пикар-3",
        "Пикар-4",
    ]

    anal = get_anal_solve(node, f_anal)
    eyler, pikar_1, pikar_2, pikar_3, pikar_4, x = get_numerical_methods(
        node, numerical, 0.5, 0
    )

    for i in range(len(anal)):
        row = [x[i], anal[i], eyler[i]]
        try:
            picar = [pikar_1[i], pikar_2[i], pikar_3[i], pikar_4[i]]
        except:
            picar = ["Ошибка", "Ошибка", "Ошибка", "Ошибка"]
        row.extend(picar)
        task.add_row(row)

    return task


def solve_third_task(node: Node) -> PrettyTable:
    task = PrettyTable()

    f_eyler = lambda x, u: np.power(x, 2) + np.power(u, 2)
    approx_picar_1 = lambda u: np.power(u, 3) / 3
    approx_picar_2 = lambda u: approx_picar_1(u) + np.power(u, 7) / 63
    approx_picar_3 = (
        lambda u: approx_picar_2(u)
        + 2 * np.power(u, 11) / 2079
        + np.power(u, 15) / 59535
    )
    approx_picar_4 = (
        lambda u: approx_picar_3(u)
        + (2 / 93555) * u**15
        + (2 / 3393495) * u**19
        + (2 / 2488563) * u**19
        + (2 / 86266215) * u**23
        + (1 / 99411543) * u**23
        + (2 / 3341878155) * u**27
        + (1 / 109876902975) * u**31
    )

    numerical = Numerical(
        f_eyler=f_eyler,
        approx_picar_1=approx_picar_1,
        approx_picar_2=approx_picar_2,
        approx_picar_3=approx_picar_3,
        approx_picar_4=approx_picar_4,
    )

    task.field_names = [
        "Аргумент",
        "Эйлер",
        "Пикар-1",
        "Пикар-2",
        "Пикар-3",
        "Пикар-4",
    ]

    eyler = solve_euler(node, 0, 0, f_eyler)
    eyler, pikar_1, pikar_2, pikar_3, pikar_4, x = get_numerical_methods(
        node, numerical, 0, 0
    )

    for i in range(len(eyler)):
        row = [x[i], eyler[i]]
        try:
            picar = [pikar_1[i], pikar_2[i], pikar_3[i], pikar_4[i]]
        except:
            picar = ["Ошибка", "Ошибка", "Ошибка", "Ошибка"]
        row.extend(picar)
        task.add_row(row)

    return task


def main():
    x_min, x_max, h = 0, 2.1, 0.01
    # x_min, x_max, h = input("Введите: минимум x, максиму x, шаг: ").split()

    node = Node(x_min=x_min, x_max=x_max, h=h)
    print("1)(u^2 + x) * u' = 1; u(1) = 0")
    print("2) 1 - 2xuu' = u^3u'; u(0.5) = 0")
    print("3) u'(x) = x^2 + u^2; u(0) = 0")
    task_num = int(input("Введите номер задания: "))
    if task_num == 1:
        task = solve_first_task(node)
        print("Решение первого задания")
    elif task_num == 2:
        task = solve_second_task(node)
        print("Решение второго задания")
    elif task_num == 3:
        task = solve_third_task(node)
        print("Решение третьего задания")

    else:
        raise Exception

    print(task)


if __name__ == "__main__":
    main()
