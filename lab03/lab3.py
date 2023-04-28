from matplotlib import pyplot as plt

LAMBDA = [
    [300, 1.36e-2],
    [500, 1.63e-2],
    [800, 1.81e-2],
    [1100, 1.98e-2],
    [2000, 2.5e-2],
    [2400, 2.74e-2]
]

K = [
    [293, 2e-2],
    [1278, 5e-2],
    [1528, 7.8e-2],
    [1677, 1e-1],
    [2000, 1.3e-1],
    [2400, 2e-1]
]

r0 = 0.35
R = 0.5
np = 1.4
sigma = 5.668e-12
T0 = 300
F0 = 100
ALPHA = 0.05

P = lambda kn_T: 4 * np ** 2 * sigma * kn_T[0] * kn_T[1] ** 3
F = lambda k: 4 * np ** 2 * sigma * T0 ** 4 * k

h = 0.3e-2  # шаг сетки
temp_eps = 1e-6
energy_eps = 1e-2

max_iters = 300


# generate grid (сетка)
def generate_grid(start, stop, step):
    grid = [start]
    x = start
    while x < stop:
        x += step
        grid.append(x)
    if grid[-1] > stop:
        grid[-1] = stop
    return grid


# linear interpolation
def generate_table_value(table_fn, point):
    lt = len(table_fn)
    i = 0
    while i < lt - 2 and table_fn[i + 1][0] < point:
        i += 1
    dX = table_fn[i + 1][0] - table_fn[i][0]
    dY = table_fn[i + 1][1] - table_fn[i][1]
    k = dY / dX
    y = table_fn[i][1]
    dx = point - table_fn[i][0]
    dy = k * dx
    return y + dy


# generate values (with table)
def generate_table_values(table_fn, points):
    return list(map(lambda x: generate_table_value(table_fn, x), points))


# generate value (with function equation)
def generate_func_values(func, points):
    return list(map(lambda x: func(x), points))


# function find f_(i-1/2)
def points_2_half_points(points):
    res = []
    for i in range(1, len(points)):
        res.append(0.5 * (points[i] + points[i - 1]))
    return res


# function Vn
def generate_vn(points, p=1):
    vn = []
    half_points = points_2_half_points(points)
    for i in range(1, len(half_points)):
        zpos = half_points[i]
        zneg = half_points[i - 1]
        v_i = (zpos ** (p + 1) - zneg ** (p + 1)) / (p + 1)
        vn.append(v_i)
    return vn


# find a, b, c, d coeefs SLAU
def generate_abcd(points, ts, p=1):
    h = points[1] - points[0]
    coef = 1 / (R ** 2 * h)
    half_points = points_2_half_points(points)
    half_kappa = points_2_half_points(generate_table_values(LAMBDA, ts)) ##
    kn = generate_table_values(K, ts)
    kn_t = list(map(lambda i: (kn[i], ts[i]), range(len(points))))
    pn = list(map(P, kn_t))
    fn = list(map(F, kn))
    vn = generate_vn(points, p)
    res = []
    a = 0
    b = half_points[0] * half_kappa[0] * coef
    b += h / 8 * (pn[0] + pn[1]) * half_points[0]
    b += h / 4 * pn[0] * points[0]
    c = half_points[0] * half_kappa[0] * coef
    c -= h / 8 * (pn[0] + pn[1]) * half_points[0]
    d = points[0] * F0 / R
    d += h / 4 * (fn[0] * points[0] + (fn[0] + fn[1]) * half_points[0] / 2)
    res.append([a, b, c, d])
    for i in range(1, len(points) - 1):
        zpos = half_points[i]
        zneg = half_points[i - 1]
        a = coef * zneg * half_kappa[i - 1]
        c = coef * zpos * half_kappa[i]
        b = a + c + pn[i] * vn[i - 1]
        d = fn[i] * vn[i - 1]
        res.append([a, b, c, d])
    a = half_points[-1] * half_kappa[-1] * coef
    a -= h / 8 * (pn[-1] + pn[-2]) * half_points[-1] * (1/R) ##
    b = points[-1] * ALPHA / R
    b += half_points[-1] * half_kappa[-1] * coef
    b += h / 8 * half_points[-1] * (pn[-1] + pn[-2]) / 2 * (1/R) ##
    b += h / 4 * pn[-1] * points[-1] * (1/R) ##
    c = 0
    d = points[-1] * ALPHA * T0 / R
    d += h / 4 * (fn[-1] + fn[-2]) * half_points[-1] * (1/R) ##
    d += h / 4 * fn[-1] * points[-1] * (1 / R) ##
    res.append([a, b, c, d])
    return res


def integrate(fn, points, ts):
    h = points[1] - points[0]
    integral = (fn(points[0], ts[0]) + fn(points[-1], ts[-1])) / 2
    for i in range(1, len(points) - 1):
        integral += fn(points[i], ts[i])
    integral *= h
    return integral


def integral_fn(z, t):
    k = generate_table_value(K, t)
    return k * (t ** 4 - T0 ** 4) * z


# метод простой итерации
def simple_iteration_method():
    points = generate_grid(r0 / R, 1, h)
    N = len(points)
    ts = [T0 for _ in points]
    p = 1
    it = 0
    while it < max_iters:
        abcd = generate_abcd(points, ts, p)
        it += 1
        ksi = [0 for _ in range(N)]
        eta = [0 for _ in range(N)]

        # прямой ход (определение прогоночных коэффициентов)
        for i in range(1, N):
            den = (abcd[i - 1][1] - abcd[i - 1][0] * ksi[i - 1])
            ksi[i] = abcd[i - 1][2] / den
            eta[i] = (abcd[i - 1][3] + abcd[i - 1][0] * eta[i - 1]) / den

        y = [0 for _ in range(N)]
        y[N - 1] = (abcd[-1][3] + abcd[-1][0] * eta[-1]) / (
                abcd[-1][1] - abcd[-1][0] * ksi[-1])
        max_yi = abs((ts[N - 1] - y[N - 1]) / ts[N - 1])

        # обратный ход (находим функцию y) - правая прогонка
        for i in range(N - 2, -1, -1):
            y[i] = y[i + 1] * ksi[i + 1] + eta[i + 1]
            max_yi = max(max_yi, abs((ts[i] - y[i]) / ts[i]))
        f1 = r0 * F0 - R * ALPHA * (y[-1] - T0)
        f2 = 4 * np ** 2 * sigma * integrate(integral_fn, points, ts) * R ** 2
        try:
            dfi = abs((f1 - f2) / f1)
            print(dfi)
        except ZeroDivisionError:
            dfi = 0
        if max_yi <= temp_eps and dfi <= energy_eps:
            break
        ts = y
    return points, ts, it, N


_points, _ts, _it, n = simple_iteration_method()
plt.figure(figsize=(6, 6))
plt.xlim(min(_points), max(_points))
plt.xlabel("r/R")
plt.ylim(0, max(_ts) * 1.5)
plt.ylabel("Температура, K")
plt.title(f"График зависимости температуры T(z)\n"
          f"F0={F0},alpha={ALPHA},h={h}, N={n}.\n"
          f"Число итераций={_it}")
plt.plot(_points, _ts)
plt.show()
