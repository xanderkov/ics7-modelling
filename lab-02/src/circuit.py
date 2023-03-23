from typing import List
import numpy as np
from src.methods import (
    find_t0_m,
    find_sigma,
    integral,

)


class Circuit:
    
    def __init__(self, R_k, H, i_arr, to_arr, m_arr, t_arr, sigma_arr):
        self.R = 0.35
        self.l_e = 12
        self.L_k = 187 * (10**(-6))
        self.C_k = 268 * (10**(-6))
        self.R_k = R_k
        self.U_co = 1400
        self.i_o = 3 #0.3
        self.T_w = 2000
        self.H = H
        self.STEP = 1e-3

        self.i_arr = i_arr
        self.t0_arr = to_arr
        self.m_arr = m_arr
        self.t_arr = t_arr
        self.sigma_arr = sigma_arr

        self.table1_filename = "./data/table1.csv"
        self.table2_filename = "./data/table2.csv"

    def di_dt(self, i: float, u: float):
        return (u - (self.R_k + self.R) * i) / self.L_k

    def t_func(self, t0: float, z: float, m: float):
        return t0 + (self.T_w - t0) * np.power(z, m)

    def r_func(self, s):
        return self.l_e / (2 * np.pi * self.R * self.R * s)

    def find_r(self, i: float):

        to, m = find_t0_m(i, self.i_arr, self.t0_arr, self.m_arr)
        t_arr2 = []
        sigma_arr2 = []
        h = self.STEP
        z = 0
        z_max = 1
        while z < z_max + h:
            t = self.t_func(to, z, m)
            sigma = find_sigma(t, self.t_arr, self.sigma_arr)
            t_arr2.append(t)
            sigma_arr2.append(sigma)
            z = z + h
        s = integral(t_arr2, sigma_arr2)
        r = self.r_func(s)
        return r

    def di_dt_func(self, i: float, u: float, r_res: float):
        return (u - (self.R_k + r_res) * i) / self.L_k

    def phi(self, i: float):
        return -(i / self.C_k)

    def runge4_with_r(self, to=0, t_max=0.01):
        i_n = self.i_o
        u_n = self.U_co
        t_n = to

        t_res = [to]
        i_res = [i_n]
        u_res = [u_n]
        r0 = self.find_r(self.i_o)
        to = find_t0_m(self.i_o, self.i_arr, self.t0_arr, self.m_arr)
        r_res = [r0]
        to_res = [to]

        while t_n < t_max:
            r_1 = self.find_r(i_n)
            k1 = self.H * self.di_dt_func(i_n, u_n, r_1)
            q1 = self.H * self.phi(i_n)

            k2 = self.H * self.di_dt_func(i_n + k1 / 2, u_n + q1 / 2, self.find_r(i_n + k1 / 2))
            q2 = self.H * self.phi(i_n + k1 / 2)

            k3 = self.H * self.di_dt_func(i_n + k2 / 2, u_n + q2 / 2, self.find_r(i_n + k1 / 2))
            q3 = self.H * self.phi(i_n + k2 / 2)

            k4 = self.H * self.di_dt_func(i_n + k3, u_n + q3, self.find_r(i_n + k3))
            q4 = self.H * self.phi(i_n + k3)

            t_n = t_n + self.H
            i_n = i_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            u_n = u_n + (q1 + 2 * q2 + 2 * q3 + q4) / 6

            r_p = self.find_r(i_n)
            to = find_t0_m(i_n, self.i_arr, self.t0_arr, self.m_arr)

            t_res.append(t_n)
            i_res.append(i_n)
            u_res.append(u_n)
            r_res.append(r_p)
            to_res.append(to)

        return t_res, i_res, u_res, r_res, to_res

    def runge4(self, t0=0, t_max=0.01):
        i_n = self.i_o
        u_n = self.U_co
        t_n = t0

        t_res = [t_n]
        i_res = [i_n]
        u_res = [u_n]

        while t_n < t_max:
            k1 = self.H * self.di_dt(i_n, u_n)
            q1 = self.H * self.phi(i_n)
            k2 = self.H * self.di_dt(i_n + k1 / 2, u_n + q1 / 2)
            q2 = self.H * self.phi(i_n + k1 / 2)
            k3 = self.H * self.di_dt(i_n + k2 / 2, u_n + q2 / 2)
            q3 = self.H * self.phi(i_n + k2 / 2)
            k4 = self.H * self.di_dt(i_n + k3, u_n + q3)
            q4 = self.H * self.phi(i_n + k3)

            t_n = t_n + self.H
            i_n = i_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            u_n = u_n + (q1 + 2 * q2 + 2 * q3 + q4) / 6

            t_res.append(t_n)
            i_res.append(i_n)
            u_res.append(u_n)

        return t_res, i_res, u_res

    def runge2_with_r(self, t0=0, t_max=0.01, beta=1 / 2):
        i_n = self.i_o
        u_n = self.U_co
        t_n = t0

        t_res = [t_n]
        i_res = [i_n]
        u_res = [u_n]
        r0 = self.find_r(i_n)
        t0 = find_t0_m(i_n, self.i_arr, self.t0_arr, self.m_arr)
        r_res = [r0]
        t0_res = [t0]

        while t_n < t_max:
            k1 = self.H * self.di_dt_func(i_n, u_n, self.find_r(i_n))
            q1 = self.H * self.phi(i_n)
            k2 = self.H * self.di_dt_func(i_n + k1 / (2 * beta), u_n + q1 / (2 * beta), self.find_r(i_n + k1 / (2 * beta)))
            q2 = self.H * self.phi(i_n + k1 / (2 * beta))

            t_n = t_n + self.H
            i_n = i_n + (1 - beta) * k1 + beta * k2
            u_n = u_n + (1 - beta) * q1 + beta * q2

            r0 = self.find_r(i_n)
            
            t0 = find_t0_m(i_n, self.i_arr, self.t0_arr, self.m_arr)

            t_res.append(t_n)
            i_res.append(i_n)
            u_res.append(u_n)
            r_res.append(r0)
            t0_res.append(t0)

        return t_res, i_res, u_res, r_res, t0_res

    def runge2(self, t0=0, t_max=0.01, beta=1 / 2):
        i_n = self.i_o
        u_n = self.U_co
        t_n = t0

        t_res = [t0]
        i_res = [i_n]
        u_res = [u_n]

        while t_n < t_max:
            k1 = self.H * self.di_dt(i_n, u_n)
            q1 = self.H * self.phi(i_n)
            k2 = self.H * self.di_dt(i_n + k1 / (2 * beta), u_n + q1 / (2 * beta))
            q2 = self.H * self.phi(i_n + k1 / (2 * beta))

            t_n = t_n + self.H
            i_n = i_n + (1 - beta) * k1 + beta * k2
            u_n = u_n + (1 - beta) * q1 + beta * q2

            t_res.append(t_n)
            i_res.append(i_n)
            u_res.append(u_n)
        return t_res, i_res, u_res

    def euler_with_r(self, t0=0, t_max=0.01):
        i_n = self.i_o
        u_n = self.U_co
        t_n = t0

        t_res = [t0]
        i_res = [i_n]
        u_res = [u_n]
        r0 = self.find_r(i_n)
        t0 = find_t0_m(i_n, self.i_arr, self.t0_arr, self.m_arr)
        r_res = [r0]
        t0_res = [t0]

        while t_n < t_max:
            k1 = self.H * self.di_dt_func(i_n, u_n, self.find_r(i_n))
            q1 = self.H * self.phi(i_n)

            t_n = t_n + self.H
            i_n = i_n + k1
            u_n = u_n + q1

            r_p = self.find_r(i_n)
            t0 = find_t0_m(i_n, self.i_arr, self.t0_arr, self.m_arr)
            t_res.append(t_n)
            i_res.append(i_n)
            u_res.append(u_n)
            r_res.append(r_p)
            t0_res.append(t0)

        return t_res, i_res, u_res, r_res, t0_res

    def euler(self, t0=0, t_max=0.01):
        i_n = self.i_o
        u_n = self.U_co
        t_n = t0

        t_res = [t0]
        i_res = [i_n]
        u_res = [u_n]

        while t_n < t_max:
            k1 = self.H * self.di_dt(i_n, u_n)
            q1 = self.H * self.phi(i_n)

            t_n = t_n + self.H
            i_n = i_n + k1
            u_n = u_n + q1

            t_res.append(t_n)
            i_res.append(i_n)
            u_res.append(u_n)

        return t_res, i_res, u_res
