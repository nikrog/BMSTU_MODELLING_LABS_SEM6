import numpy as np
import numpy as np
import matplotlib.pyplot as plt
from calc import solve
from parameters import a, b, hx, hz, N, M


def plot_graph_3d(z, title):
    fig = plt.figure(figsize=(12, 10))
#    plt.title(title)
    ax = fig.add_subplot(111, projection='3d')
    a_x = 0
    b_x = b
    h_x = hz
    a_y = 0
    b_y = a
    h_y = hx
    x = np.linspace(a_x, b_x, N)
    y = np.linspace(a_y, b_y, M)
    x, y = np.meshgrid(x, y)
    ax.plot_surface(x, y, z)
    ax.set_xlabel('Z')
    ax.set_ylabel('X')
    ax.set_zlabel('u(X, Z)')
    ax.set_zlim(200, 3000)
    plt.show()


def plot_graphs_2d(u, fixed_x=1, step=1):
    plt.figure(figsize=(6, 6))
    if fixed_x:
        plt.xlabel("z, см")
        plt.ylabel("T, K")
        plt.title(f"Функции u(x0, z)\n")
        a_z = 0
        b_z = b
        z_arr = np.linspace(a_z, b_z, M)
        for x in [0, N//2, N-1]: #range(0, len(u), step):
            plt.plot(z_arr, u[x], label=f"u({x*hx:.3f},z)")
    else:
        plt.xlabel("x, см")
        plt.ylabel("T, K")
        plt.title(f"Функции u(x, z0)\n")
        a_x = 0
        b_x = a
        x_arr = np.linspace(a_x, b_x, N)
        for z in [0, M//2, M-1]: #range(0, M//2, step):
            c = [r[z] for r in u]
            plt.plot(x_arr, c,  label=f"u(x,{z*hz:.3f})")
    plt.legend()
    plt.show()


def main():
    one_only = input('one only (y/n):').strip() != 'n'
    x = np.linspace(0, a, N)
    y = np.linspace(0, b, M)
    X, Y = np.meshgrid(x, y)
    t, solution = solve(one_only)
    z = np.array(solution)
    Z = z.reshape(X.shape)
    plot_graph_3d(Z, "Решение задачи")
    plot_graphs_2d(Z, 0, 5)
    plot_graphs_2d(Z, 1, 30)


if __name__ == '__main__':
    main()
