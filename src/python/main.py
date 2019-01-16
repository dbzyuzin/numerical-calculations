from mymetods import *
from task import *
import numpy as np


maxiter = 100
maxiter_jacobi = 333

ys = np.zeros((N1, N2))
Cs = [[list() for _ in range(N2)] for _ in range(N1)]
F = np.zeros((N1, N2))

edge_computing(ys, u, x1, x2)

for iter_count in range(maxiter):
    #iter coef
    for i in range(1, N1-1):
        for j in range(1, N2-1):
            Csi = np.zeros(5)
            Csi[0] = (h2/h1*(ki(ys[i+1, j], ys[i, j]) + ki(ys[i-1, j], ys[i, j])) + \
                      + h1/h2*(ki(ys[i, j+1], ys[i, j]) + \
                      + ki(ys[i, j-1], ys[i, j])) +\
                      + h1*h2*q(ys[i,j]))
            Csi[1] = h2/h1*ki(ys[i+1, j], ys[i, j])
            Csi[2] = h2/h1*ki(ys[i-1, j], ys[i, j])
            Csi[3] = h1/h2*ki(ys[i, j+1], ys[i, j])
            Csi[4] = h1/h2*ki(ys[i, j-1], ys[i, j])

            F[i,j] = h1*h2*f(ys[i, j])
            Cs[i][j] = Csi

    if test_solution(ys, Cs, F) < eps: break

    #jacobi
    for iter_count_j in range(maxiter_jacobi):
        for i in range(1, N1-1):
            for j in range(1, N2-1):
                ys[i, j] = (F[i,j] + Cs[i][j][1]*ys[i+1, j] + Cs[i][j][2]*ys[i-1, j] + \
                    + Cs[i][j][3]*ys[i, j+1] + Cs[i][j][4]*ys[i, j-1])/Cs[i][j][0]
        if test_solution(ys, Cs, F) < eps_j:
            break

#Задаем точное решение
ysol = solution(u, x1, x2)

#вычисляем ошибку на сетке
yerr = final_error(ys, ysol)

# выводим количество итераций и максимальную ошибку
print(make_res(N1, N2, h1, h2, eps, iter_count, yerr))
with open("test.txt", "a") as f:
    f.write(make_res(N1, N2, h1, h2, eps, iter_count, yerr, False))
