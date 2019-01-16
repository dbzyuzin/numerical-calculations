from mymetods import *
from task import *
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation



maxiter = 100
maxiter_jacobi = 20

ys = np.zeros((N1, N2))
Cs = [[list() for _ in range(N1)] for _ in range(N2)]
F = np.zeros((N1, N2))

edge_computing(ys, u, x1, x2)

for iter_count in range(maxiter):
    #iter coef
    for i in range(1,N1-1):
        for j in range(1,N2-1):
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

    #testing
    rka = 0.0
    for i in range(1, N1-1):
        for j in range(1, N2-1):
            rka = max(rka, abs(F[i,j] + Cs[i][j][1]*ys[i+1, j] + \
                      + Cs[i][j][2]*ys[i-1, j] + Cs[i][j][3]*ys[i, j+1] + \
                      + Cs[i][j][4]*ys[i, j-1] - Cs[i][j][0]*ys[i, j]))

    if rka < eps: break

    #jacobi

    for iter_count_j in range(maxiter_jacobi):
        for i in range(1, N1-1):
            for j in range(1, N2-1):
                ys[i, j] = (F[i,j] + Cs[i][j][1]*ys[i+1, j] + Cs[i][j][2]*ys[i-1, j] + Cs[i][j][3]*ys[i, j+1] + Cs[i][j][4]*ys[i, j-1])/Cs[i][j][0]

#Задаем точное решение
ysol = np.zeros((N1, N2))
yerr = ys.copy()
for i in range(N1):
    for j in range(N2):
        ysol[i,j] = u(x1[i],x2[j])
#вычисляем ошибку на сетке
rka = 0.0
for i in range(N1):
    for j in range(N2):
        yerr[i, j] = abs(ysol[i, j]-ys[i, j])

# выводим количество итераций и максимальную ошибку
print(make_res(N1, N2, h1, h2, eps, iter_count, yerr))
with open("test.txt", "a") as f:
    f.write(make_res(N1, N2, h1, h2, eps, iter_count, yerr, False))
#заготовка построения графиков
def my3d_plot(zdata, x1, x2):
    X, Y = np.meshgrid(x1, x2)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, zdata, color='0.75', rstride=1, cstride=1, cmap="magma")
    ax.view_init(elev=10., azim=0)

my3d_plot(ys, x1, x2)
my3d_plot(ysol, x1, x2)
