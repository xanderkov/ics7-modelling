import numpy as np
# T ( z) = T0 + (Tw - T0 ) z^m
def find_dx(x1, x2):
    return x1 - x2


def find_t0_m(i, i_arr, to_arr, m_arr):
    n = len(i_arr)
    j = 0

    if i < i_arr[0]:
        m = m_arr[0]
        to = to_arr[0]
        return to, m
    elif i > i_arr[n - 1]:
        m = m_arr[n - 1]
        to = to_arr[n - 1]
        return to, m

    while i_arr[j] > i or j == n - 2:
        j += 1
    j -= 1

    if j < n - 1:
        dx = i_arr[j + 1] - i_arr[j]
        di = i - i_arr[j]
        to = to_arr[j] + ((to_arr[j + 1] - to_arr[j]) * di / dx)
        m = m_arr[j] + ((m_arr[j + 1] - m_arr[j]) * di / dx)
        # print(j, i, i_arr[j+1], i_arr[j])
    else:
        dx = i_arr[n - 1] - i_arr[n - 2]
        di = i - i_arr[n - 1]
        to = to_arr[n - 2] + ((to_arr[n - 1] - to_arr[n - 2]) * di / dx)
        m = m_arr[n - 1]

    if m < 0:
        print(i_arr[-1])
        # print(m, i, fl)
    return to, m


def find_sigma(t, t_arr, sigma_arr):
    n = len(t_arr)
    j = 0
    if t < t_arr[0]:
        sigma = sigma_arr[0]
        return sigma

    elif t > t_arr[n - 1]:
        sigma = sigma_arr[n]
        return sigma

    while True:
        if t_arr[j] > t or j == n - 2:
            break
        j += 1
    j -= 1

    # while j < n - 1 and t_arr[j] > t:
    #   j += 1

    if j < n - 1:
        dx = t_arr[j+1] - t_arr[j]
        di = t - t_arr[j]
        sigma = sigma_arr[j] + ((sigma_arr[j + 1] - sigma_arr[j]) * di / dx)
        # print(t, t_arr[j+1], t_arr[j])
    else:
        dx = find_dx(t_arr[n - 1], t_arr[n - 2])
        di = t - t_arr[n - 1]
        sigma = sigma_arr[n - 2] + ((sigma_arr[n - 1] - sigma_arr[n - 2]) * di / dx)

    return sigma


# Trapetions
def integral(arr1, arr2):
    l = len(arr1)
    s = 0
    for i in range(l - 1):
        s += ((arr2[i] + arr2[i + 1]) / 2) * (arr1[i + 1] - arr1[i])
    return s


def interpolation(x_arr, y_arr, h):
    x_arr = [np.log(x) for x in x_arr]
    y_arr = [np.log(y) for y in y_arr]
    res_x = []
    res_y = []
    for i in range(len(x_arr) - 1):
        dx = x_arr[i + 1] - x_arr[i]
        dy = y_arr[i + 1] - y_arr[i]
        k = dy / dx
        x = x_arr[i]
        y = y_arr[i]
        # k *= y / x
        while (x + h) < x_arr[i + 1]:
            res_x.append(x)
            res_y.append(y)
            x += h
            y += k * h
        if i == len(x_arr) - 1:
            res_x.append(x)
            res_y.append(y)
    res_x = [np.exp(x) for x in res_x]
    res_y = [np.exp(y) for y in res_y]
    return res_x, res_y
