import numpy as np

def make_res(N1, N2, h1, h2, eps, iter_count, yerr, flag=True):
    if flag:
        s = "Параметры :\n{: ^7}{: ^7}{: ^7}{: ^7}{: ^8}\n".format("N1", "N2", "h1", "h2", "eps")
        s += "{: ^7d}{: ^7d}{: ^7.3f}{: ^7.3f}{: ^8}\n\n".format(N1, N2, h1, h2, eps)
        s += "Результаты:\n {: ^12}{: ^24}\n".format("Iter count", "Max Fail")
        s += " {: ^12}{: ^24.15f}\n{}\n".format(iter_count, np.max(yerr)), "="*36)
    else:
        s = "{: ^7}{: ^7}{: ^7}{: ^7}{: ^8}".format("N1", "N2", "h1", "h2", "eps")
        s += "{: ^12}{: ^24}\n".format("Iter count", "Max Fail")
        s += "{: ^7d}{: ^7d}{: ^7.3f}{: ^7.3f}{: ^8}".format(N1, N2, h1, h2, eps)
        s += " {: ^12}{: ^24.15f}\n{}\n".format(iter_count, np.max(yerr)), "="*72)
    return s

def edge_computing(ys, u, x1, x2):
    N1, N2 = np.shape(ys)
    for i, x in enumerate(x1):
        ys[i, 0] = u(x,0)
        ys[i, -1] = u(x,1)
    for i, x in enumerate(x2):
        ys[0, i] = u(0, x)
        ys[-1, i] = u(1, x)

def final_error(theory, solution):
    eps = 0
    N, M = np.shape(theory)
    for i in range(0, N):
        for j in range(0, M):
            d = np.abs(solution[i,j]-theory[i,j])
            if d > eps:
                eps = d
    return eps
