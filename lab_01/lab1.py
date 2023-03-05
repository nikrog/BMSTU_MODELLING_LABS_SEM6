from prettytable import PrettyTable
from math import e

STEP_TABLE = 0.01
MAX_X = 2.01
STEP = 1e-5


def f1(u, x):
    return u * u + x


def f1p1(u):
    return 1 + u + pow(u, 3) / 3


def f1p2(u):
    return 1 + u + pow(u, 2) / 2 + pow(u, 3) / 3 + pow(u, 4) / 12


def f1p3(u):
    return 1 + u + pow(u, 2) / 2 + pow(u, 3) / 2 + pow(u, 4) / 12 + pow(u, 5) / 60


def f1p4(u):
    return 1 + u + pow(u, 2) / 2 + pow(u, 3) / 2 + pow(u, 4) / 8 + pow(u, 5) / 60 + pow(u, 6) / 360


def f1_analyt(u):
    return 3 * pow(e, u) - u * u - 2 * u - 2


def f2(u, x):
    return pow(u, 3) + 2 * x * u


def f2p1(u):
    return 1 / 2 + pow(u, 2) / 2 + pow(u, 4) / 4


def f2p2(u):
    return 1 / 2 + pow(u, 2) / 2 + pow(u, 4) / 2 + pow(u, 6) / 12


def f2p3(u):
    return 1 / 2 + pow(u, 2) / 2 + pow(u, 4) / 2 + pow(u, 6) / 6 + pow(u, 8) / 48


def f2p4(u):
    return 1 / 2 + pow(u, 2) / 2 + pow(u, 4) / 2 + pow(u, 6) / 6 + pow(u, 8) / 24 + pow(u, 10) / 240


def f2_analyt(u):
    return pow(e, u * u) - (u * u) / 2 - 1 / 2


def f3(x, u):
    return x * x + u * u


def f3p1(x):
    return pow(x, 3) / 3


def f3p2(x):
    return pow(x, 3) / 3 + pow(x, 7) / 63


def f3p3(x):
    return pow(x, 3) / 3 + pow(x, 7) / 63 + (2 * pow(x, 11)) / 2079 + \
        pow(x, 15) / 59535


def f3p4(x):
    return pow(x, 3) / 3 + pow(x, 7) / 63 + (2 * pow(x, 11)) / 2079 + \
        (13 * pow(x, 15)) / 218295 + (82 * pow(x, 19)) / 37328445 + \
        (662 * pow(x, 23)) / (23 * 453835305) + \
        (4 * pow(x, 27)) / (27 * 123773265) + \
        pow(x, 31) / (31 * 59535 * 59535)


def euler(x_max, h, func, x0, y0):
    result = list()
    k = 0
    #print(y0)
    while x0 < x_max:
        #print(y0)
        if abs(x0 - k) < 1e-4:
            result.append(round(y0, 6))
            k = k + STEP_TABLE
        y0 = y0 + h * func(x0, y0)
        x0 += h
    return result


def picard(x_max, h, func, x, y):
    result = list()
    k = 0
    while x < x_max:
        if abs(x - k) < 1e-4:
            result.append(round(y, 6))
            k = k + STEP_TABLE
        x += h
        y = func(x)
    return result


def analytical(x_max, h, func, x, y):
    result = list()
    k = 0
    while x < x_max:
        if abs(x - k) < 1e-4:
            result.append(round(y, 6))
            k = k + STEP_TABLE
        x += h
        y = func(x)
    return result


def x_range(x_max, h):
    result = list()
    x = 0
    k = 0
    while x < x_max:
        if abs(x - k) < 1e-4:
            result.append(round(x, 2))
            k = k + STEP_TABLE
        x += h
    return result


def main():
    tb = PrettyTable()
    task = int(input("Input task number: "))
    if task == 1:
        tb.add_column("x", x_range(MAX_X, STEP))
        tb.add_column("Аналит.", analytical(MAX_X, STEP, f1_analyt, 0, 1))
        tb.add_column("Эйлера (явный)", euler(MAX_X, STEP, f1, 0, 1))
        tb.add_column("Пикара 1", picard(MAX_X, STEP, f1p1, 0, 1))
        tb.add_column("Пикара 2", picard(MAX_X, STEP, f1p2, 0, 1))
        tb.add_column("Пикара 3", picard(MAX_X, STEP, f1p3, 0, 1))
        tb.add_column("Пикара 4", picard(MAX_X, STEP, f1p4, 0, 1))
    elif task == 2:
        tb.add_column("x", x_range(MAX_X, STEP))
        tb.add_column("Аналит.", analytical(MAX_X, STEP, f2_analyt, 0, 0.5))
        tb.add_column("Эйлера (явный)", euler(MAX_X, STEP, f2, 0, 0.5))
        tb.add_column("Пикара 1", picard(MAX_X, STEP, f2p1, 0, 0.5))
        tb.add_column("Пикара 2", picard(MAX_X, STEP, f2p2, 0, 0.5))
        tb.add_column("Пикара 3", picard(MAX_X, STEP, f2p3, 0, 0.5))
        tb.add_column("Пикара 4", picard(MAX_X, STEP, f2p4, 0, 0.5))
    elif task == 3:
        tb.add_column("x", x_range(MAX_X, STEP))
        tb.add_column("Эйлера (явный)", euler(MAX_X, STEP, f3, 0, 0))
        tb.add_column("Пикара 1", picard(MAX_X, STEP, f3p1, 0, 0))
        tb.add_column("Пикара 2", picard(MAX_X, STEP, f3p2, 0, 0))
        tb.add_column("Пикара 3", picard(MAX_X, STEP, f3p3, 0, 0))
        tb.add_column("Пикара 4", picard(MAX_X, STEP, f3p4, 0, 0))
    else:
        print("Incorrect task number!")
        return
    print(tb)



if __name__ == "__main__":
    main()
