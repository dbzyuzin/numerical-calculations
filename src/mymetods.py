from numpy import array, zeros, diag, diagflat, dot

def jacobi(A, f, maxIter=225, x=None, eps=1e-3):
    """
    A - заданная матрица, f - столбец справа от нуля,
    eps -  предпологаемая погрешность, x - начальное приблежение.
    """
    if x is None:
        x = zeros(len(A[0]))

    D = diag(A)
    R = A - diagflat(D)

    for i in range(maxIter):
        x = (f - dot(R,x)) / D
        if all(abs(dot(A, x) - f) < eps):
            break

    return x


def jacobi_test():
    eps = 1e-4
    res = True

    A = array([[10.0,1.0, -1],[1.0,10.0, -1], [-1, 1, 10]], dtype=float)
    f = array([11.0,10.0, 10], dtype=float)
    sol = jacobi(A,f, eps=eps)
    if any(abs(dot(A, sol) - f) >= eps): res = False


    A = array([[2.0,1.0],[5.0,7.0]], dtype=float)
    f = array([11.0,13.0])
    sol = jacobi(A,f, eps=eps)
    if any(abs(dot(A, sol) - f) >= eps): res = False

    print("Test jacobi() Success")
    return res
