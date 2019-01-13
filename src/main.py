from mymetods import *
from task import *
import numpy as np

N = 100
x1 = np.linspace(0, 1, N)
x2 = np.linspace(0, 1, N)
h = 1/N
eps = 10**-4
ys = np.ones((N, N))*0.99

ytest = solution(x1, x2)
ys = np.zeros((N, N))

edge_computing(ys, u, (x1, x2))
A = first_appr(ys, (N,N))
tau, mu = simple_iter_koef(A)
ns = simple_iter_count(eps, mu)
A, F = make_koef_matr(ys, (N, N), tau, (h,h))

X = np.zeros((N*N,1))
X0 = np.zeros((N*N,1))
z = 0
for i in range(0, N):
    for j in range(0, N):
        X[z,0] = ys[i,j]
        X0[z,0]= ytest[i,j]
        z = z+1
while all(abs(X0-X)) > eps:
    X = jacobi(A,F, ns, X, eps)
    for z in range(0, N*N):
        ys[z // N, z % N] = X[z]
    tau, mu = simple_iter_koef(A)
    ns = simple_iter_count(eps, mu)
    A, F = make_koef_matr(ys, (N,N), tau, (h,h))
