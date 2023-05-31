from math import exp

# Physics params
a = 10
b = 10
N = 100
M = 100
u0 = 300
t_min = 0
t_step = 0.01
t_n = 100
a1 = 0.0134
b1 = 1
c1 = 4.35 * 10 ** -4
m1 = 1
f0 = 30#10
beta = 1.5#0.1
alpha2 = 10#0.05
alpha3 = 10#0.05
alpha4 = 10#0.05
F0 = 0#30

# User params
x0 = a / 2
z0 = b / 2


# Derived params
def l(u):
    return a1 * (b1 + c1 * u ** m1)


hx = a / (N - 1)
hz = b / (M - 1)


def f(x_idx, z_idx):
    x = x_idx * hx
    z = z_idx * hz
    return f0 * exp(-beta * (x - x0) ** 2 * (z - z0) ** 2)


# Introduced params
def c(u):
    return 1


K = [
    [293, 2e-2],
    [1278, 5e-2],
    [1528, 7.8e-2],
    [1677, 1e-1],
    [2000, 1.3e-1],
    [2400, 2e-1]
]

# Tech params
max_x_iters = 200
max_z_iters = 200
eps = 1e-4
