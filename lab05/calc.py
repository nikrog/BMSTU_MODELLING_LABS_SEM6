import parameters as prm


def generate_table_value(table_fn, point):
    # lt = len(table_fn)
    # i = 0
    # while i < lt - 2 and table_fn[i + 1][0] < point:
    #     i += 1
    # dX = table_fn[i + 1][0] - table_fn[i][0]
    # dY = table_fn[i + 1][1] - table_fn[i][1]
    # k = dY / dX
    # b = table_fn[i][1]
    # dx = point - table_fn[i][0]
    # dy = k * dx
    # return b + dy
    return 0.05


def do_sweep(abcd, N):
    ksi = [0 for _ in range(N)]
    eta = [prm.u0 for _ in range(N)]

    for i in range(2, N):
        den = (abcd[i - 1][1] - abcd[i - 1][0] * ksi[i - 1])
        ksi[i] = abcd[i - 1][2] / den
        eta[i] = (abcd[i - 1][3] + abcd[i - 1][0] * eta[i - 1]) / den

    y = [prm.u0 for _ in range(N)]
    #y[N - 1] = (abcd[-1][3] + abcd[-1][0] * eta[-1]) / (
            #abcd[-1][1] - abcd[-1][0] * ksi[-1])
    for i in range(N - 2, -1, -1):
        y[i] = y[i + 1] * ksi[i + 1] + eta[i + 1]
    return y


def solve_for_t_at_z(prev_u_matrix, _t, m):
    k_left = generate_table_value(prm.K, prev_u_matrix[0][m])
    k_right = generate_table_value(prm.K, prev_u_matrix[-1][m])
    abcd = [
        #[0, -k_left / prm.hx, -k_left / prm.hx, -prm.F0]  # x = 0
        [0, -1, 0, -prm.u0]
    ]
    # first 0 is not used
    l = [0] + [prm.l(prev_u_matrix[n][m]) for n in range(1, prm.N)]
    for n in range(1, prm.N - 1):
        c = prm.c(prev_u_matrix[n][m])
        K1 = l[n] / prm.hx ** 2
        K3 = l[n + 1] / prm.hx ** 2
        K2 = K1 + K3 + 2 * c / prm.t_step
        K4 = 2 * c * prev_u_matrix[n][m] / prm.t_step + prm.f(n, m)
        if 0 < m < prm.M - 1:
            K4 += (l[n + 1] - l[n]) * (
                    prev_u_matrix[n][m + 1] - prev_u_matrix[n][m]
            ) / prm.hz ** 2
            K4 += l[n] * (prev_u_matrix[n][m + 1] - 2 * prev_u_matrix[n][m] +
                          prev_u_matrix[n][m - 1]) / prm.hz ** 2
        abcd.append([K1, K2, K3, K4])
    abcd.append(
        #[k_right / prm.hx, k_right / prm.hx, 0, prm.alpha2 * (prm.u0 - prev_u_matrix[-1][m])]  # x = a
        [0, -1, 0, prm.u0]
    )
    u_column = do_sweep(abcd, prm.N)
    check = check_convergence_list(
        [prev_u_matrix[n][m] for n in range(prm.N)],
        u_column
    )
    cur_u_matrix = [
        [
            prev_u_matrix[i][j] for j in range(prm.M)
        ]
        for i in range(prm.N)
    ]
    for n in range(prm.N):
        cur_u_matrix[n][m] = u_column[n]
    return check, cur_u_matrix


def solve_for_t_at_x(prev_u_matrix, _t, n):
    k_left = generate_table_value(prm.K, prev_u_matrix[n][0])
    k_right = generate_table_value(prm.K, prev_u_matrix[n][-1])
    abcd = [
        # [
        #     0,
        #     k_left / prm.hz,
        #     k_left / prm.hz,
        #     prm.alpha3 * (prm.u0 - prev_u_matrix[n][0])
        # ]  # z = 0
        [0, -1, 0, prm.u0]
    ]
    # first 0 is not used
    l = [0] + [prm.l(prev_u_matrix[n][m]) for m in range(1, prm.M)]
    for m in range(1, prm.M - 1):
        c = prm.c(prev_u_matrix[n][m])
        K1 = l[m] / prm.hz ** 2
        K3 = l[m + 1] / prm.hz ** 2
        K2 = K1 + K3 + 2 * c / prm.t_step
        K4 = 2 * c * prev_u_matrix[n][m] / prm.t_step + prm.f(n, m)
        if 0 < n < prm.N - 1:
            K4 += (l[m + 1] - l[m]) * (
                    prev_u_matrix[n + 1][m] - prev_u_matrix[n][m]
            ) / prm.hz ** 2
            K4 += l[m] * (prev_u_matrix[n + 1][m] - 2 * prev_u_matrix[n][m] +
                          prev_u_matrix[n - 1][m]) / prm.hz ** 2
        abcd.append([K1, K2, K3, K4])
    abcd.append(
        # [
        #     k_right / prm.hz,
        #     k_right / prm.hz,
        #     0,
        #     prm.alpha4 * (prm.u0 - prev_u_matrix[n][-1])
        # ]  # z = b
        [0, -1, 0, prm.u0]
    )
    u_row = do_sweep(abcd, prm.M)
    check = check_convergence_list(prev_u_matrix[n], u_row)
    cur_u_matrix = [
        [
            prev_u_matrix[i][j] for j in range(prm.M)
        ]
        for i in range(prm.N)
    ]
    for m in range(prm.M):
        cur_u_matrix[n][m] = u_row[m]
    return check, cur_u_matrix


def solve_for_t(prev_u_matrix, t, var, file=None):
    cur_u_matrix = prev_u_matrix
    if var == "x":
        for m in range(prm.M):
            it = 0
            while it < prm.max_x_iters:
                it += 1
                check, cur_u_matrix = solve_for_t_at_z(cur_u_matrix, t, m)
                if check:
                    break
            if m % 10 == 0:
                print(f"done m = {m} for t = {t}, it = {it}")
    else:
        for n in range(prm.N):
            it = 0
            while it < prm.max_z_iters:
                it += 1
                check, cur_u_matrix = solve_for_t_at_x(cur_u_matrix, t, n)
                if check:
                    break
            if n % 10 == 0:
                print(f"done n = {n} for t = {t}, it = {it}")
    return cur_u_matrix


def check_convergence_list(prev_list, next_list):
    max_diff = 0
    for i in range(len(prev_list)):
        cur_diff = abs(prev_list[i] - next_list[i]) / abs(prev_list[i])
        if cur_diff > prm.eps:
            return False
        max_diff = max(
            max_diff,
            cur_diff
        )
    return max_diff <= prm.eps


def check_convergence(prev_u_matrix, next_u_matrix, file=None):
    max_diff = 0
    for n in range(prm.N):
        for m in range(prm.N):
            cur_diff = abs(prev_u_matrix[n][m] - next_u_matrix[n][m]) / abs(
                prev_u_matrix[n][m])
            if cur_diff > prm.eps:
                print(cur_diff, file=file)
                return False
            max_diff = max(
                max_diff,
                cur_diff
            )
    if max_diff <= prm.eps:
        print(max_diff, file=file)
    return max_diff <= prm.eps


def solve(only_one=False):
    u_matrix = [[prm.u0 for _ in range(prm.M)] for _ in range(prm.N)]
    t = prm.t_min
    file = open("f.txt", "w")
    for i in range(prm.t_n):
        u_matrix_half = solve_for_t(u_matrix, t + prm.t_step / 2, "x", file)
        u_matrix_next = solve_for_t(u_matrix_half, t + prm.t_step, "z", file)
        if check_convergence(u_matrix, u_matrix_next, file):
            break
        u_matrix = [row[:] for row in u_matrix_next]
        t += prm.t_step
        if only_one:
            break
    file.close()
    return t, u_matrix


if __name__ == "__main__":
    print(list(solve(True)))
