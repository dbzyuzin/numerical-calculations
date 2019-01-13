import numpy as np

def k(u):
    return 1.0 + u**2

def ki(u1, u2):
    k1 = k(u1)
    k2 = k(u2)
    return 2*(k1*k2)/(k1+k2)

def q(u):
    return 0.5*np.pi**2*(1-u**2)

def u(x1, x2):
    return np.cos(0.5*np.pi*(x1+x2))

def f(u):
    return np.pi**2*u**3

def solution(x1, x2):
    sol = np.ones((len(x1), len(x2)))

    for i in range(0, len(x1)):
        for j in range(0, len(x2)):
            sol[i,j] = u(x1[i],x2[j])
    return sol

def make_koef_matr(ys, size, tau, h):
    N, M = size
    h1, h2 = h
    Asize = N*M
    A = np.zeros((Asize, Asize))
    F = np.zeros((Asize, 1))
    for i in range(0, Asize):
        ys1 = ys[i,i]
        ki1 = ki(ys[i+1,i], ys1)
        ki2 = ki(ys[i-1,i], ys1)
        kj1 = ki(ys[i,i+1], ys1)
        kj2 = ki(ys[i,i-1], ys1)
        km1 = tau*(ki1 + ki2 + kj1 + kj2)**-1
        A[i,i] = 1 - tau
        A[i+1,i] = ki1*km1
        A[i-1,i] = ki2*km1
        A[i,i+1] = kj1*km1
        A[i,i-1] = kj2*km1
        F[i,1] = (h1*h2 * f(ys1) - q(ys1)*ys1)*km1
    return (A, F)

def first_appr(ys, size):
    N, M = size
    K = np.ones((N+1,N+1))
    for i in range(0, N):
        for j in range(0, M):
            ys1 = ys[i,j]
            ki1 = ki(ys[i+1,j], ys1)
            ki2 = ki(ys[i-1,j], ys1)
            kj1 = ki(ys[i,j+1], ys1)
            kj2 = ki(ys[i,j-1], ys1)
            K[i,j] = (ki1 + ki2 + kj1 + kj2)
    return K