import numpy as np

N1, N2 = 20, 20
x1 = np.linspace(0, 1, N1)
x2 = np.linspace(0, 1, N2)
h1 = 1/N1
h2 = 1/N2
eps = 1e-5


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
