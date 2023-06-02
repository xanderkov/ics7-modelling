import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import MaxNLocator

f0 = 30
betta = 5
x0 = 5
z0 = 5
T0 = 300
dx = 1e-1
dz = 1e-1

x1 = 0
z1 = 0
x2 = 10
z2 = 10

a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1

tau = 0.5

def f(x: float, z: float) -> float:
    return f0 * math.exp(-betta * (math.pow((x - x0), 2) * math.pow((z - z0), 2))) # + f0 * math.exp(-betta * (math.pow((x - 2), 2) + math.pow((z - 2), 2)))

def lambda_func(u: float) -> float:
    # return 1.0
    return a1 * (b1 + c1 * math.pow(u, m1))

def Ax() -> float:
    return dz * dz * tau / 2

def Bx() -> float:
    return dx * dx * dz * dz + tau * dz * dz

def Cx() -> float:
    return dz * dz * tau / 2

def Dx(u_1: float, u: float, u1: float, x: float, z: float) -> float:
    return dx * dx * dz * dz * u + tau / 2 * dx * dx * (u_1  - 2 * u + u1) + tau * dx * dx * dz * dz * f(x, z) / 2 / lambda_func(u)
    # -f(x, z) * dx * dx * dz * dz / lambda_func(u) - dx * dx * (u_1 + u1)

def Az() -> float:
    return dx * dx * tau / 2

def Bz() -> float:
    return dx * dx * dz * dz + tau * dx * dx

def Cz() -> float:
    return dx * dx * tau / 2

def Dz(u_1: float, u: float, u1: float, x: float, z: float) -> float:
    return dx * dx * dz * dz * u + tau / 2 * dz * dz * (u_1  - 2 * u + u1) + tau * dx * dx * dz * dz * f(x, z) / 2 / lambda_func(u)

def border() -> float:
    return T0

def interpolate(xval: float, xs: list[float], ys: list[float]) -> float:
    i = 0
    while (i < len(xs) and xval < xs[i]):
        i += 1
    if (i >= len(xs) - 1):
        i = len(xs) - 2
    yval = ys[i] + (ys[i + 1] - ys[i])  / (xs[i + 1] - xs[i]) * (xval - xs[i])
    return yval

def method(a: list[float], b: list[float], c: list[float], d: list[float]) -> list[float]:
    y: list[float] = [0] * len(b)
    xi: list[float] = [0] * (len(b) - 1)
    eta: list[float] = [0] * (len(b) - 1)
    xi[0] = 0
    eta[0] = T0
    for i in range(1, len(xi)):
        tmp = b[i] - a[i] * xi[i - 1]
        xi[i] = c[i] / tmp
        eta[i] = (a[i] * eta[i - 1] + d[i]) / tmp

    y[len(b) - 1] = T0
    for i in range(0, len(y) - 1):
        num = len(y) - 1 - i
        y[num - 1] = xi[num - 1] * y[num] + eta[num - 1]

    return y

def form_coefficients_x(z: float, grid: list[list[float]]) -> tuple[list[float], list[float], list[float], list[float]]:
    x = x1
    x_index = 0
    z_index = int((z - z1) / dz)
    As: list[float] = []
    Bs: list[float] = []
    Cs: list[float] = []
    Ds: list[float] = []
    while x < x2 - 1e-5:
        As.append(Ax())
        Bs.append(Bx())
        Cs.append(Cx())
        Ds.append(Dx(grid[x_index][z_index - 1], grid[x_index][z_index], grid[x_index][z_index + 1], x, z))
        x += dx
        x_index += 1

    return (As, Bs, Cs, Ds)

def form_coefficients_z(x: float, grid: list[list[float]]) -> tuple[list[float], list[float], list[float], list[float]]:
    z = z1
    z_index = 0
    x_index = int((x - x1) / dx)
    As: list[float] = []
    Bs: list[float] = []
    Cs: list[float] = []
    Ds: list[float] = []
    while z < z2 - 1e-5:
        As.append(Az())
        Bs.append(Bz())
        Cs.append(Cz())
        Ds.append(Dz(grid[x_index - 1][z_index], grid[x_index][z_index], grid[x_index + 1][z_index], x, z))
        z += dz
        z_index += 1

    return (As, Bs, Cs, Ds)

def iteration(grid: list[list[float]]) -> tuple[list[list[float]], float]:
    nx = len(grid) - 1
    nz = len(grid[0]) - 1
    new_grid = []
    for row in grid:
        new_grid.append(row.copy())
    max_err = 5000

    for i in range(0, nx + 1):
        grid[i][0] = T0
        grid[i][nz] = T0
    
    for i in range(0, nz + 1):
        grid[0][i] = T0
        grid[nx][i] = T0
    
    its = 0
    while (its < 1 and max_err > 1e-2):
        max_err = 0
        its += 1
        for i in range(1, nz):
            z = z1 + i * dz
            As, Bs, Cs, Ds = form_coefficients_x(z, grid)
            us = method(As, Bs, Cs, Ds)
            for k in range(0, nx + 1):
                new_grid[k][i] = us[k]
                max_err = max(max_err, math.fabs(new_grid[k][i] - grid[k][i]))

    its = 0
    while (its < 1 and max_err > 1e-2):
        max_err = 0
        its += 1
        for i in range(1, nx):
            x = x1 + i * dx
            As, Bs, Cs, Ds = form_coefficients_z(x, grid)
            us = method(As, Bs, Cs, Ds)
            for k in range(0, nz + 1):
                new_grid[i][k] = us[k]
                max_err = max(max_err, math.fabs(new_grid[i][k] - grid[i][k]))

    return (new_grid, max_err)

def run() -> list[list[float]]:
    err = 5e100
    its = 0
    row: list[float] = [T0] * int((z2 - z1) / dz)
    grid: list[list[float]] = [row.copy()] * int((x2 - x1) / dx)
    while its < 1 and err >= 1e-2:
        (grid, err) = iteration(grid)
        its += 1
        print(err)

    return grid

grid = run()

X = []
x = x1
while x < x2 - 1e-5:
    X.append(x)
    x += dx

Z = []
z = z1
while z < z2 - 1e-5:
    Z.append(z)
    z += dz

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

y = tuple(grid)
xgrid, zgrid = np.meshgrid(X, Z)
ygrid = np.array(grid, np.float64)
ax.plot_surface(xgrid, zgrid, ygrid, color='r', alpha=0.4)
ax.set_xlabel('X axis')
ax.set_zlabel('Y axis')
ax.set_ylabel('Z axis')
plt.show()
