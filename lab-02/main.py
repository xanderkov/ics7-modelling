from src.circuit import Circuit
import pandas as pd
from prettytable import PrettyTable
import matplotlib.pyplot as plt


def menu():
    print("1) Calculate by Euler")
    print("2) Calculate by Runge2")
    print("3) Calculate by Runge4")


def main():
    table1 = pd.read_csv("./data/table1.csv")
    table2 = pd.read_csv("./data/table2.csv")
    answer = PrettyTable()

    H = 1e-4
    MIN = 0.01
    MAX = 600 * 1e-6

    i_arr = table1["I"]
    t0_arr = table1["T0"]
    m_arr = table1["m"]
    t_arr = table2["T"]
    sigma_arr = table2["Delta"]

    circuit = Circuit(i_arr, t0_arr, m_arr, t_arr, sigma_arr)

    menu()

    choose = int(input("Choose numbers (1 - 3): "))
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=False, sharey=False, figsize=(10, 10))
    ax1, ax2, ax3, ax4 = axes.flatten()
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

    if choose == 1:
        x_res, y_res, z_res, r_res, t_res = circuit.euler_with_r()

        ax1.set_title("I(t) RK1")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK1")
        ax2.plot(x_res, z_res, colors[1])
        ax3.set_title("Rp(t) RK1")
        ax3.plot(x_res, r_res, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        ax4.set_title("I(t) * Rp(t) RK1")
        ax4.plot(x_res, ir_res, colors[3])

    elif choose == 2:
        x_res, y_res, z_res, r_res, t_res = circuit.runge2_with_r()

        ax1.set_title("I(t) RK2")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK2")
        ax2.plot(x_res, z_res, colors[1])
        ax3.set_title("Rp(t) RK2")
        ax3.plot(x_res, r_res, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        ax4.set_title("I(t) * Rp(t) RK2")
        ax4.plot(x_res, ir_res, colors[3])

    elif choose == 3:
        x_res, y_res, z_res, r_res, t_res = circuit.runge4_with_r()

        ax1.set_title("I(t) RK4")
        ax1.plot(x_res, y_res, colors[0])
        ax2.set_title("U(t) RK4")
        ax2.plot(x_res, z_res, colors[1])
        ax3.set_title("Rp(t) RK4")
        ax3.plot(x_res, r_res, colors[2])
        ir_res = [y_res[i] * r_res[i] for i in range(len(y_res))]
        ax4.set_title("I(t) * Rp(t) RK4")
        ax4.plot(x_res, ir_res, colors[3])

    plt.show()

    

if __name__ == "__main__":
    main()
    