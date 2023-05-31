import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

M = 4
eps = 10 ** (-3)

# функция расчета коэффициентов для системы dxdt
def FindCofficientODE(x_centre, param):
    x_m = np.zeros((2, M + 1))
    for i in range(0, 2):
        x_m[i][0] = x_centre[i]
    for j in range(1, M + 1):
        for i in range(0, 2):
            x_m[i][j] += x_m[i][j-1] * param[i][0]
            summa1 = 0
            for p in range(0, j):
                summa1 += x_m[0][p] * x_m[1][j - 1 - p]
            x_m[i][j] -= param[i][1] * summa1
            summa2 = 0
            for p in range(0, j):
                summa2 += x_m[i][p] * x_m[i][j - 1 - p]
            x_m[i][j] -= param[i][2] * summa2
            x_m[i][j] /= j
    return x_m

# функция расчета коэффициентов для системы dydt
def FindCofficientVar(y_centre, param, x12):
    y_m = np.zeros((12, M + 1))
    for i in range(0, 12):
        y_m[i][0] = y_centre[i]
    a1 = param[0][0] - param[0][1] *x12[1] - 2 * param[0][2] * x12[0]
    a2 = -1 * param[0][1] * x12[0]
    a3 = param[1][0] - param[1][1] * x12[0] - 2 * param[1][2] * x12[1]
    a4 = -1 * param[1][1] * x12[1]
    a_m = np.array([[a1, a2], [a3, a4]])
    for j in range(1, M + 1):
        for i in range(0, 6):
            y_m[i][j] = (a_m[0][0] * y_m[i][j - 1] + a_m[0][1] * y_m[i + 6][j - 1])/(j + 1)
        for i in range(6, 12):
            y_m[i][j] = (a_m[1][0] * y_m[i][j - j] + a_m[1][1] * y_m[i - 6][j - 1])/(j + 1)
    y_m[0][1] += (x12[0]) / (j + 1)
    y_m[1][1] += (-1 * x12[0] * x12[1]) / (j + 1)
    y_m[2][1] += (-1 * x12[0] * x12[0]) / (j + 1)
    y_m[9][1] += (x12[1]) / (j + 1)
    y_m[10][1] += (-1 * x12[0] * x12[1]) / (j + 1)
    y_m[11][1] += (-1 * x12[1] * x12[1]) / (j + 1)
    return y_m

# функция определения шага
def ChooseStepODE(x_centre, hn, par):
    x0_max = x_centre.max()
    Q = np.zeros(2)
    for i in range(0, 1):
        Q[i] = abs(par[i][0] * x0_max + par[i][1] * x0_max ** (2) + par[i][2] * x0_max ** (2))
    Qmax = Q.max()
    rx = 1 / Qmax
    h_i = max(rx * 0.9, hn)
    d = abs(h_i) / rx
    error = x0_max * (1 - d) ** (-1) * (d) ** (M + 1)
    while error > 10 ** (-6):
        h_i /= 2
        d = abs(h_i) / rx
        error = x0_max * (1 - d) ** (-1) * (d) ** (M + 1)
    return h_i

def ChooseStepVar(h, par, y, x):
    s1 = abs(par[0][0] - par[0][1] * x[1] - 2 * par[0][2] * x[0]) + abs(par[0][1] * x[0])
    s2 = abs(par[1][0] - par[1][1] * x[0] - 2 * par[1][2] * x[1]) + abs(par[1][1] * x[1])
    s = max(s1, s2)
    ry = 1/s
    a = math.sqrt(pow(x[0], 2) + pow(x[0], 4) + pow(-x[0]*x[1], 2) + pow(x[1], 2) + pow(x[1], 4))
    y0_abs = max(abs(y))
    h = min(h, 0.9 * ry)
    d = abs(h) / ry
    err = (y0_abs + a * ry) * math.exp(d) * (d ** (M + 1)) / math.factorial(M + 1)
    while err > 10 ** (-3):
        h /= 2.0
        d = abs(h) / ry
        err = ((y0_abs + a * ry) * math.exp(d) * (d ** (M + 1))) / math.factorial(M + 1)
    return h

# функция расчета ряда Тейлора
def TaylorSeriesODE(xm, step):
    taylor_x = np.zeros(2)
    for i in range(0, 2):
        for m in range (0, M + 1):
            taylor_x[i] += (xm[i][m] * step ** (m)) / math.factorial(m)
    return taylor_x

def TaylorSeriesVar(ym, step):
    taylor_y = np.zeros(12)
    for i in range(0, 12):
        for m in range (0, M + 1):
            taylor_y[i] += (ym[i][m] * step ** (m)) / math.factorial(m)
    return taylor_y

def TaylorLV(par, x_centre):
    x = np.zeros((2, T))
    x[0][0] = x_centre[0]
    x[1][0] = x_centre[1]
    coef = np.zeros((2, M + 1))
    taylor_series = np.zeros(2)
    x0 = np.array([x_centre[0], x_centre[1]])
    h = 0
    for t in range(1, T):
        l = 0.03
        while (h < t):
            l = ChooseStepODE(x0, l, par)
            h += l
            if h > t:
                l = t - h
                h = t
            coef = FindCofficientODE(x0, par)
            taylor_series = TaylorSeriesODE(coef, l)
            x0[0] = taylor_series[0]
            x0[1] = taylor_series[1]
        x[0][t] = x0[0]
        x[1][t] = x0[1]
    return x

def TaylorVariation(par, y_centre, x12):
    coef = np.zeros((12, M + 1))
    taylor_series = np.zeros(12)
    y = np.zeros((12, T))
    y0 = y_centre
    h = 0.0
    for t in range(1, T):
        l = 0.1
        while (h < t):
            l = ChooseStepVar(l, par, y0, x12[:, t - 1])
            h += l
            if h > t:
                l = h - t
                h = t
            coef = FindCofficientVar(y0, par, x12[:, t - 1])
            taylor_series = TaylorSeriesVar(coef, l)
            for j in range(0, 12):
                y0[j] = taylor_series[j]
        for i in range(0, 12):
            y[i][t] = y0[i]
    return y


def calculate_Taylor(k_start, data_x, years):
    # решаем систему Лотки-Вольтерры с начальным приближением
    data_x0 = np.array([data_x[0][0], data_x[1][0]])
    global T
    T = years.size
    x_sol = TaylorLV(k_start, data_x0)
    functional = 0
    for i in range(0, 2):
        for j in range(0, T):
            functional += 1/2.0 * ((x_sol[i][j] - data_x[i][j]) ** (2))
    fun = np.zeros(2)
    fun[0] = functional

    # решаем систему уравнений в вариациях
    y0 = np.zeros(12)
    y_sol = TaylorVariation(k_start, y0, x_sol)

    # рассчитываем градиентные уравнения
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


    x0 = np.array([data_x0[0], data_x0[1]])
    functional = 0
    x_sol = TaylorLV(k_new, x0)
    for i in range(0, 2):
        for j in range(0, T):
            functional += 1/2.0 * ((x_sol[i][j] - data_x[i][j]) ** (2))
    fun[1] = functional   

    while (abs(fun[1]- fun[0]) > s2):
        y0 = np.zeros(12)
        y_sol = TaylorVariation(k_new, y0, x_sol)

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
        x_sol = TaylorLV(k_new, x0)
        for i in range(0, 2):
            for j in range(0, T):
                functional += 1/2.0 * ((x_sol[i][j] - data_x[i][j]) ** (2))
        fun[0] = fun[1]
        fun[1] = functional
        if (fun[1] >= fun[0]):
            s2 /= 2
    x_sol = TaylorLV(k_new, data_x0)
    # plt.plot(years, x_sol[0], years, x_sol[1])
    # plt.show()
    return x_sol[0], x_sol[1],  k_new






