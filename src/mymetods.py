import numpy as np

def jacobi(A, f, maxIter=225, x=None, eps=1e-3):
    """
    A - заданная матрица, f - столбец справа от нуля,
    eps -  предпологаемая погрешность, x - начальное приблежение.
    """
    if x is None:
        x = np.zeros(len(A[0]))

    D = np.diag(A)
    R = A - np.diagflat(D)

    for _ in range(maxIter):
        x = (f - np.dot(R,x)) / D
        if all(abs(np.dot(A, x) - f) < eps):
            break

    return x


def jacobi_test():
    eps = 1e-4
    res = True

    A = np.array([[10.0,1.0, -1],[1.0,10.0, -1], [-1, 1, 10]], dtype=float)
    f = np.array([11.0,10.0, 10], dtype=float)
    sol = jacobi(A,f, eps=eps)
    if any(abs(np.dot(A, sol) - f) >= eps): res = False


    A = np.array([[2.0,1.0],[5.0,7.0]], dtype=float)
    f = np.array([11.0,13.0])
    sol = jacobi(A,f, eps=eps)
    if any(abs(np.dot(A, sol) - f) >= eps): res = False

    print("Test jacobi() Success")
    return res

def simple_iter_koef(A):
    eigs, _ = np.linalg.eig(np.abs(A))
    lambmin = np.min(eigs)
    lambmax = np.max(eigs)
    return ( 2/(lambmax+lambmin) , lambmax/lambmin )

def simple_iter_count(eps, mu):
    return int(1 + np.round(-np.log(eps)/np.log((mu+1)/(mu-1))))

def edge_computing(ys, u, X):
    x1, x2 = X
    N, M = np.shape(ys)
    for i in range(0, N):
        ys[i, 0] = u(x1[i],0)
        ys[i, -1] = u(x1[i],1)
    for i in range(0, M):
        ys[0, i] = u(0, x2[i])
        ys[-1, i] = u(1, x2[i])

def final_error(theory, solution):
    eps = 0
    N, M = np.shape(theory)
    for i in range(0, N):
        for j in range(0, M):
            d = np.abs(solution[i,j]-theory[i,j])
            if d > eps:
                eps = d
    return eps