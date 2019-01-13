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

jacobi(A,F, ns, None, eps)
