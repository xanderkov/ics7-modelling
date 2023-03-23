from src.circuit import Circuit
import pandas as pd
from prettytable import PrettyTable
import matplotlib.pyplot as plt


def menu():
    print("1) Calculate by Euler")
    print("2) Calculate by Runge2")
    print("3) Calculate by Runge4")
    print("4) Calculate by Runge4 Rk=200")
    print("5) Calculate by Runge4 Rk+Rp=0")


def main():
    table1 = pd.read_csv("./data/table1.csv")
    table2 = pd.read_csv("./data/table2.csv")
    answer = PrettyTable()

    MIN = 0.01
    MAX = 600 * 1e-6

    i_arr = list(table1.I.array)
    t0_arr = list(table1.T0.array)
    m_arr = list(table1.m.array)
    t_arr = list(table2.TK.array)
    sigma_arr = list(table2.Delta.array)

    circuit = Circuit(0.25, i_arr, t0_arr, m_arr, t_arr, sigma_arr)

    menu()

    choose = int(input("Choose numbers (1 - 5): "))
    fig, axes = plt.subplots(nrows=2, ncols=3, sharex=False, sharey=False, figsize=(10, 10))
    ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:black']

    if choose == 1:
        x_res, y_res, z_res, r_res, t_res = circuit.euler_with_r(0, MAX)

        ax1.set_title("I(t) RK1")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK1")
        ax2.plot(x_res, z_res, colors[1])
        ax3.set_title("Rp(t) RK1")
        ax3.plot(x_res, r_res, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        ax4.set_title("I(t) * Rp(t) RK1")
        ax4.plot(x_res, ir_res, colors[3])
        ax5.set_title("T0(t) RK1")
        ax5.plot(x_res, t_res, colors[4])

    elif choose == 2:
        x_res, y_res, z_res, r_res, t_res = circuit.runge2_with_r(0, MAX)

        ax1.set_title("I(t) RK2")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK2")
        ax2.plot(x_res, z_res, colors[1])
        ax3.set_title("Rp(t) RK2")
        ax3.plot(x_res, r_res, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        ax4.set_title("I(t) * Rp(t) RK2")
        ax4.plot(x_res, ir_res, colors[3])
        ax5.set_title("T0(t) RK2")
        ax5.plot(x_res, t_res, colors[4])

    elif choose == 3:
        x_res, y_res, z_res, r_res, t_res = circuit.runge4_with_r(0, MAX)

        ax1.set_title("I(t) RK4")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK4")
        ax2.plot(x_res, z_res, colors[1])
        ax3.set_title("Rp(t) RK4")
        ax3.plot(x_res, r_res, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        ax4.set_title("I(t) * Rp(t) RK4")
        ax4.plot(x_res, ir_res, colors[3])
        ax5.set_title("T0(t) RK4")
        ax5.plot(x_res, t_res, colors[4])

    elif choose == 4:
        circuit = Circuit(200, i_arr, t0_arr, m_arr, t_arr, sigma_arr)
        x_res, y_res, z_res, r_res, t_res = circuit.runge4_with_r(0, 20 * 1e-6)
        plt.title("I(t) (R_k = 200) RK4")
        plt.plot(x_res, y_res, colors[0])

    elif choose == 5:
        x_res, y_res, z_res = circuit.runge4(0, MAX)
        plt.title("I(t) RK4 Rp+Rk=0")
        plt.plot(x_res, y_res, colors[0])

    plt.show()

    

if __name__ == "__main__":
    main()
    