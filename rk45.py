import numpy as np

eps = 10 ** (-3)

def system_LK(x, k):
    dx1_dt = k[0][0] * x[0] - k[0][1] * x[0] * x[1] - k[0][2] * x[0] * x[0]
    dx2_dt = k[1][0] * x[1] - k[1][1] * x[1] * x[0] - k[1][2] * x[1] * x[1]
    return dx1_dt, dx2_dt

def variation_LK(y, x, k):
    dy11dt = (k[0][0] - k[0][1] * x[1] - 2 * k[0][2] * x[0]) * y[0] - k[0][1] * x[0] * y[6] + x[0]
    dy12dt = (k[0][0] - k[0][1] * x[1] - 2 * k[0][2] * x[0]) * y[1] - k[0][1] * x[0] * y[7] - x[0] * x[1]
    dy13dt = (k[0][0] - k[0][1] * x[1] - 2 * k[0][2] * x[0]) * y[2] - k[0][1] * x[0] * y[8] - x[0] * x[0]
    dy14dt = (k[0][0] - k[0][1] * x[1] - 2 * k[0][2] * x[0]) * y[3] - k[0][1] * x[0] * y[9]
    dy15dt = (k[0][0] - k[0][1] * x[1] - 2 * k[0][2] * x[0]) * y[4] - k[0][1] * x[0] * y[10]
    dy16dt = (k[0][0] - k[0][1] * x[1] - 2 * k[0][2] * x[0]) * y[5] - k[0][1] * x[0] * y[11]
    dy21dt = (k[1][0] - k[1][1] * x[0] - 2 * k[1][2] * x[1]) * y[6] - k[1][1] * x[1] * y[0]
    dy22dt = (k[1][0] - k[1][1] * x[0] - 2 * k[1][2] * x[1]) * y[7] - k[1][1] * x[1] * y[1]
    dy23dt = (k[1][0] - k[1][1] * x[0] - 2 * k[1][2] * x[1]) * y[8] - k[1][1] * x[1] * y[2]
    dy24dt = (k[1][0] - k[1][1] * x[0] - 2 * k[1][2] * x[1]) * y[9] - k[1][1] * x[1] * y[3] + x[1]
    dy25dt = (k[1][0] - k[1][1] * x[0] - 2 * k[1][2] * x[1]) * y[10] - k[1][1] * x[1] * y[4] - x[0] * x[1]
    dy26dt = (k[1][0] - k[1][1] * x[0] - 2 * k[1][2] * x[1]) * y[11] - k[1][1] * x[1] * y[5] - x[1] * x[1]
    return dy11dt, dy12dt, dy13dt, dy14dt, dy15dt, dy16dt, dy21dt, dy22dt, dy23dt, dy24dt, dy25dt, dy26dt

# def find_coefficient(x, h, param):
#     q = np.zeros((2, 6))
#     for i in range(0, 2):
#         x0 = x
#         q[i][0] = h * system_LK(x0, param)[i]
#         x0 = x + 1.0/4.0 * q[i][0]
#         q[i][1] = h * system_LK(x0, param)[i]
#         x0 = x + 3.0/32.0 * q[i][0] + 9.0/32.0 * q[i][1]
#         q[i][2] = h * system_LK(x0, param)[i]
#         x0 = x + 1932.0/2197.0 * q[i][0] - 7200.0/2197.0 * q[i][1] + 7296.0/2197.0 * q[i][2]
#         q[i][3] = h * system_LK(x0, param)[i]
#         x0 = x + 439.0/216.0 * q[i][0] - 8.0 * q[i][1] + 3680.0/513.0 * q[i][2] - 845.0/4104.0 * q[i][3]
#         q[i][4] = h * system_LK(x0, param)[i]
#         x0 = x - 8.0/27.0 * q[i][0] + 2.0 * q[i][1] - 3544.0/2565.0 * q[i][2] + 1859.0/4104.0 * q[i][3] - 11.0/40.0 * q[i][4]
#         q[i][5] = h * system_LK(x0, param)[i]
#     return q

# def step_check(q):
#     E1 = 1.0/360.0 * q[0][0] - 128.0/4275.0 * q[0][2] - 2197.0/75240.0 * q[0][3] + 1.0/50.0 * q[0][4] + 2.0/55.0 * q[0][5]
#     E2 = 1.0/360.0 * q[1][0] - 128.0/4275.0 * q[1][2] - 2197.0/75240.0 * q[1][3] + 1.0/50.0 * q[1][4] + 2.0/55.0 * q[1][5]
#     return max(abs(E1), abs(E2))

def step_check_RK4(q, m):
    E = np.zeros(m)
    for i in range(0, m):
        E[i] = 1/6 * (q[i][0] - 4 * q[i][1] + 2 * q[i][2] + q[i][3])
    return np.max(abs(E))

def find_RK4_ode(x, h, param):
    q = np.zeros((2, 4))
    for i in range(0, 2):
        x0 = x
        q[i][0] = h * system_LK(x0, param)[i]
        x0 = x + 0.5 * q[i][0]
        q[i][1] = h * system_LK(x0, param)[i]
        x0 = x + 0.5 * q[i][1]
        q[i][2] = h * system_LK(x0, param)[i]
        x0 = x + q[i][2]
        q[i][3] = h * system_LK(x0, param)[i]
    return q

def find_RK4_var(y, h, param, xv):
    q = np.zeros((12, 4))
    for i in range(0, 12):
        y0 = y
        q[i][0] = h * variation_LK(y0, xv, param)[i]
        y0 = y + 0.5 * q[i][0]
        q[i][1] = h * variation_LK(y0, xv, param)[i]
        y0 = y + 0.5 * q[i][1]
        q[i][2] = h * variation_LK(y0, xv, param)[i]
        y0 = y + q[i][2]
        q[i][3] = h * variation_LK(y0, xv, param)[i]
    return q

def RK_ode(param, x_0):
    coef = np.zeros((2, 4))
    len = 0.0
    x = np.zeros((2, T))
    for i in range(0, 2):
        x[i][0] = x_0[i]
    for t in range(1, T):
        h = 0.03
        while (len < t):
            coef = find_RK4_ode(x_0, h, param)
            E = step_check_RK4(coef, 2)
            while (E > eps):
                h /= 2
                coef = find_RK4_ode(x_0, h, param)
                E = step_check_RK4(coef, 2)
            if (len + h) > t:
                h = t - len
                coef = find_RK4_ode(x_0, h, param)
                E = step_check_RK4(coef, 2)
            len += h
            if (E < eps/8.0):
                h *= 2
            for i in range(0, 2):
                x_0[i] += 1/6 * (coef[i][0] + 2 * coef[i][1] + 2 * coef[i][2] + coef[i][3])
        for i in range(0, 2):
            x[i][t] = x_0[i]
    return x

def RK_var(param, y_0, xf):
    coef = np.zeros((12, 4))
    h = 0.0
    len = 0.0
    y = np.zeros((12, T))
    for i in range(0, 12):
        y[i][0] = y_0[i]
    for t in range(1, T):
        h = 0.1
        while (len < t):
            coef = find_RK4_var(y_0, h, param, xf[:, t - 1])
            E = step_check_RK4(coef, 12)
            while (E > eps):
                h /= 2
                coef = find_RK4_var(y_0, h, param, xf[:, t - 1])
                E = step_check_RK4(coef, 12)
            if (len + h) > t:
                h = t - len
                coef = find_RK4_var(y_0, h, param, xf[:, t - 1])
                E = step_check_RK4(coef, 12)
            len += h
            if (E < eps/8.0):
                h *= 2
            for i in range(0, 12):
                y_0[i] += 1/6 * (coef[i][0] + 2 * coef[i][1] + 2 * coef[i][2] + coef[i][3])
        for i in range(0, 12):
            y[i][t] = y_0[i]
    return y


def calculate_RK(k_start, data_x, years):
    data_x0 = np.array([data_x[0][0], data_x[1][0]])
    global T
    T = years.size
    functional = 0
    x0 = np.array([data_x0[0], data_x0[1]])
    x_sol = RK_ode(k_start, x0)
    fun = np.zeros(2)

    for i in range(0, 2):
        for j in range(0, T):
            functional += 1/2.0 * ((x_sol[i][j] - data_x[i][j]) ** (2))
    fun[0] = functional
    # plt.plot(years, x_sol[0], years, x_sol[1])
    # plt.show()

    y0 = np.zeros(12)
    y_sol = RK_var(k_start, y0, x_sol)

    gradient = np.zeros(6)
    for p in range(0, 6):
        for i in [0]:
            for j in range(0, T):
                gradient[p] -= (x_sol[i][j] - data_x[i][j]) * y_sol[i + p][j]
        for i in [1]:
            for j in range(0, T):
                gradient[p] -= (x_sol[i][j] - data_x[i][j]) * y_sol[i + p + 5][j]
    k_find = np.zeros((6, 3))
    k_find[0][0] = k_start[0][0]
    k_find[1][0] = k_start[0][1]
    k_find[2][0] = k_start[0][2]
    k_find[3][0] = k_start[1][0]
    k_find[4][0] = k_start[1][1]
    k_find[5][0] = k_start[1][2]
    s1 = 0.0000001
    s2 = 10 ** (-12)

    for i in range(0, 6):
        k_find[i][1] = k_find[i][0] + s1 * gradient[i]

    k_new = np.array([[k_find[0][1], k_find[1][1], k_find[2][1]], [k_find[3][1], k_find[4][1], k_find[5][1]]])


    functional = 0
    x0 = np.array([data_x0[0], data_x0[1]])
    x_sol = RK_ode(k_new, x0)
    for i in range(0, 2):
        for j in range(0, T):
            functional += 1/2.0 * ((x_sol[i][j] - data_x[i][j]) ** (2))
    fun[1] = functional   

    while (abs(fun[1]-fun[0]) > s2):
        y0 = np.zeros(12)
        y_sol = RK_var(k_new, y0, x_sol)

        gradient = np.zeros(6)
        for p in range(0, 6):
            for i in [0]:
                for l in range(0, T):
                    gradient[p] -= (x_sol[i][l] - data_x[i][l]) * y_sol[i + p][l]
            for i in [1]:
                for l in range(0, T):
                    gradient[p] -= (x_sol[i][l] - data_x[i][l]) * y_sol[i + p + 5][l]

        for n in range(0, 6):
            k_find[n][2] = k_find[n][1] + s1 * gradient[n] + 0.5 * (k_find[n][1] - k_find[n][0])
        for n in range(0, 6):
            k_find[n][0] = k_find[n][1]
        for n in range(0, 6):
            k_find[n][1] = k_find[n][2]
        k_new = np.array([[k_find[0][2], k_find[1][2], k_find[2][2]], [k_find[3][2], k_find[4][2], k_find[5][2]]])

   
        functional = 0
        x0 = np.array([data_x0[0], data_x0[1]])
        x_sol = RK_ode(k_new, x0)
        for i in range(0, 2):
            for j in range(0, T):
                functional += 1/2.0 * ((x_sol[i][j] - data_x[i][j]) ** (2))
        fun[0] = fun[1]
        fun[1] = functional
        if (fun[1] >= fun[0]):
            s2 /= 2
    x_sol = RK_ode(k_new, data_x0)
    # plt.plot(years, x_sol[0], years, x_sol[1])
    # plt.show()
    return x_sol[0], x_sol[1], k_new
