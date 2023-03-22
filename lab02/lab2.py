from math import log, exp, pi
from prettytable import PrettyTable
from matplotlib import pyplot

R = 0.35
l = 12
L_k = 187 * (10**(-6))
C_k = 268 * (10**(-6))
R_k = 0.25
R_k2 = - 0.35
U_co = 1400
I_o = 3
T_w = 2000
STEP = 1e-4


def I_analyt():
    pass


def f(t, I, U):
    return (U - (R_k + R) * I) / L_k


def f2(t, I, U, I_arr, T0_arr, m_arr, T_arr, sigma_arr):
    return (U - (R_k + find_R(I, I_arr, T0_arr, m_arr, T_arr, sigma_arr)) * I) / L_k


def phi(t, I):
    return -(I / C_k)


def runge4_withR(I0, U0, h, I_arr, T0_arr, m_arr, T_arr, sigma_arr, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f2(t_n, I_n, U_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q1 = h * phi(t_n, I_n)
        k2 = h * f2(t_n + h / 2, I_n + k1 / 2, U_n + q1 / 2, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q2 = h * phi(t_n + h / 2, I_n + k1 / 2)
        k3 = h * f2(t_n + h / 2, I_n + k2 / 2, U_n + q2 / 2, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q3 = h * phi(t_n + h / 2, I_n + k2 / 2)
        k4 = h * f2(t_n + h, I_n + k3, U_n + q3, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q4 = h * phi(t_n + h / 2, I_n + k3)

        t_n = t_n + h
        I_n = I_n + (k1 + 2*k2 + 2*k3 + k4) / 6
        U_n = U_n + (q1 + 2*q2 + 2*q3 + q4) / 6

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


def runge4(I0, U0, h, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f(t_n, I_n, U_n)
        q1 = h * phi(t_n, I_n)
        k2 = h * f(t_n + h / 2, I_n + k1 / 2, U_n + q1 / 2)
        q2 = h * phi(t_n + h / 2, I_n + k1 / 2)
        k3 = h * f(t_n + h / 2, I_n + k2 / 2, U_n + q2 / 2)
        q3 = h * phi(t_n + h / 2, I_n + k2 / 2)
        k4 = h * f(t_n + h, I_n + k3, U_n + q3)
        q4 = h * phi(t_n + h / 2, I_n + k3)

        t_n = t_n + h
        I_n = I_n + (k1 + 2*k2 + 2*k3 + k4) / 6
        U_n = U_n + (q1 + 2*q2 + 2*q3 + q4) / 6

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


def runge2_withR(I0, U0, h, I_arr, T0_arr, m_arr, T_arr, sigma_arr, t0=0, t_max=0.01, beta=1/2):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f2(t_n, I_n, U_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q1 = h * phi(t_n, I_n)
        k2 = h * f2(t_n + h / (2*beta), I_n + k1 / (2*beta), U_n + q1 / (2*beta),
                    I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q2 = h * phi(t_n + h / (2*beta), I_n + k1 / (2*beta))

        t_n = t_n + h
        I_n = I_n + (1 - beta) * k1 + beta * k2
        U_n = U_n + (1 - beta) * q1 + beta * q2

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


def runge2(I0, U0, h, t0=0, t_max=0.01, beta=1/2):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f(t_n, I_n, U_n)
        q1 = h * phi(t_n, I_n)
        k2 = h * f(t_n + h / (2*beta), I_n + k1 / (2*beta), U_n + q1 / (2*beta))
        q2 = h * phi(t_n + h / (2*beta), I_n + k1 / (2*beta))

        t_n = t_n + h
        I_n = I_n + (1 - beta) * k1 + beta * k2
        U_n = U_n + (1 - beta) * q1 + beta * q2

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


def euler_withR(I0, U0, h, I_arr, T0_arr, m_arr, T_arr, sigma_arr, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f2(t_n, I_n, U_n, I_arr, T0_arr, m_arr, T_arr, sigma_arr)
        q1 = h * phi(t_n, I_n)

        t_n = t_n + h
        I_n = I_n + k1
        U_n = U_n + q1

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


def euler(I0, U0, h, t0=0, t_max=0.01):
    I_n = I0
    U_n = U0
    t_n = t0

    t_res = [t0]
    I_res = [I0]
    U_res = [U0]

    while t_n < t_max:
        k1 = h * f(t_n, I_n, U_n)
        q1 = h * phi(t_n, I_n)

        t_n = t_n + h
        I_n = I_n + k1
        U_n = U_n + q1

        t_res.append(t_n)
        I_res.append(I_n)
        U_res.append(U_n)

    return t_res, I_res, U_res


# функция T
def T_func(T0, z, m):
    return T0 + (T_w - T0) * z**m

# функция R
def R_func(S):
    return l/(2*pi*R**2 * S)


def find_T0_m(I, I_arr, T0_arr, m_arr):
    dx = I_arr[1] - I_arr[0]
    fl = 0
    n = len(I_arr)
    j = 0
    while j < n - 1 and I_arr[j] > I:
        j += 1
    if j < n - 1:
        di = I - I_arr[j]
        T0 = T0_arr[j] + ((T0_arr[j + 1] - T0_arr[j]) * di / dx)
        m = m_arr[j] + ((m_arr[j + 1] - m_arr[j]) * di / dx)
        fl = 1
    '''
    for i in range(n - 1):
        if I > I_arr[i]:
            di = I - I_arr[i]
            T0 = T0_arr[i] + ((T0_arr[i+1] - T0_arr[i]) * di / dx)
            m = m_arr[i] + ((m_arr[i + 1] - m_arr[i]) * di / dx)
            fl = i
    '''
    if fl == 0:
        di = I - I_arr[n - 1]
        T0 = T0_arr[n - 2] + ((T0_arr[n - 1] - T0_arr[n - 2]) * di / dx)
        #m = m_arr[n - 2] + ((m_arr[n - 1] - m_arr[n - 2]) * di / dx)
        m = m_arr[n - 1]

    if m < 0:
        print(I_arr[-1])
        print(m, I, fl)
    return T0, m


def find_sigma(T, T_arr, sigma_arr):
    dx = T_arr[1] - T_arr[0]
    fl = 0
    n = len(T_arr)
    j = 0
    while j < n - 1 and T_arr[j] > T:
        j += 1
    if j < n - 1:
        di = T - T_arr[j]
        sigma = sigma_arr[j] + ((sigma_arr[j + 1] - sigma_arr[j]) * di / dx)
        fl = 1
    '''
    for i in range(n - 1):
        if T > T_arr[i]:
            di = T - T_arr[i]
            sigma = sigma_arr[i] + ((sigma_arr[i + 1] - sigma_arr[i]) * di / dx)
            fl = 1
    '''
    if fl == 0:
        di = T - T_arr[n - 1]
        sigma = sigma_arr[n - 2] + ((sigma_arr[n - 1] - sigma_arr[n - 2]) * di / dx)

    return sigma


def find_R(I, I_arr, T0_arr, m_arr, T_arr, sigma_arr):
    T0, m = find_T0_m(I, I_arr, T0_arr, m_arr)
    T_arr2 = []
    sigma_arr2 = []
    h = STEP
    z = 0
    z_max = 1
    while z < z_max + h:
        T = T_func(T0, z, m)
        sigma = find_sigma(T, T_arr, sigma_arr)
        T_arr2.append(T)
        sigma_arr2.append(sigma)
        z = z + h
    S = integral(T_arr2, sigma_arr2)
    R = R_func(S)
    return R


# метод трапеций
def integral(arr1, arr2):
    l = len(arr1)
    s = 0
    for i in range(l - 1):
        s += ((arr2[i] + arr2[i + 1]) / 2) * (arr1[i + 1] - arr1[i])
    return s


# линейная с вырвнивающими коэффициентами lnx, lny
def interpolation(x_arr, y_arr, h):
    x_arr = [log(x) for x in x_arr]
    y_arr = [log(y) for y in y_arr]
    res_x = []
    res_y = []
    for i in range(len(x_arr) - 1):
        dx = x_arr[i + 1] - x_arr[i]
        dy = y_arr[i + 1] - y_arr[i]
        k = dy / dx
        x = x_arr[i]
        y = y_arr[i]
        #k *= y / x
        while (x + h) < x_arr[i + 1]:
            res_x.append(x)
            res_y.append(y)
            x += h
            y += k * h
        if i == len(x_arr) - 1:
            res_x.append(x)
            res_y.append(y)
    res_x = [exp(x) for x in res_x]
    res_y = [exp(y) for y in res_y]
    return res_x, res_y


if __name__ == "__main__":
    tb = PrettyTable()
    I = [0.5, 1, 5, 10, 50, 200, 400, 800, 1200]
    T0 = [6730, 6790, 7150, 7270, 8010, 9185, 10010, 11140, 12010]
    m = [0.50, 0.55, 1.7, 3, 11, 32, 40, 41, 39]
    T = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000]
    sigma = [0.031, 0.27, 2.05, 6.06, 12.0, 19.9, 29.6, 41.1, 54.1, 67.7, 81.5]

    #res_x, res_y = interpolation(I, m, 0.0001)
    iI, iT0 = interpolation(I, T0, 0.0001)
    iI, im = interpolation(I, m, 0.0001)
    print(len(iI), len(iT0), len(im))
    iT, isigma = interpolation(T, sigma, 0.0001)
    x_res, y_res, z_res = euler_withR(I_o, U_co, 0.0001, I, T0, m, T, sigma, 0, 0.01)
    x_res2, y_res2, z_res2 = runge2_withR(I_o, U_co, 0.0001, I, T0, m, T, sigma, 0, 0.01)
    x_res3, y_res3, z_res3 = runge4_withR(I_o, U_co, 0.0001, I, T0, m, T, sigma, 0, 0.01)
    x_res = [round(x, 4) for x in x_res]
    y_res = [round(x, 6) for x in y_res]
    y_res2 = [round(x, 6) for x in y_res2]
    y_res3 = [round(x, 6) for x in y_res3]
    tb.add_column("x", x_res)
    tb.add_column("Эйлера", y_res)
    tb.add_column("Рунге 2", y_res2)
    tb.add_column("Рунге 4", y_res3)
    print(tb)
    #x_res, y_res, z_res = euler(I_o, U_co, 0.0001, 0, 0.01)
    #x_res, y_res, z_res = runge2(I_o, U_co, 0.0001, 0, 0.01)
    #x_res, y_res, z_res = runge4(I_o, U_co, 0.0001, 0, 0.02)
    #print("x:", x_res)
    #print("y:", z_res)

    print(x_res)
    #pyplot.plot(iI, im)
    pyplot.plot(x_res, y_res)
    pyplot.plot(x_res, z_res)
    #pyplot.plot(x_res, z_res)
    #pyplot.plot(x_res, z_res)
    #pyplot.plot(iT, isigma)
    pyplot.show()