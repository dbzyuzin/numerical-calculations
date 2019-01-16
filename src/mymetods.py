import numpy as np

def make_res(N1, N2, h1, h2, eps, iter_count, yerr, flag=True):
    if flag:
        s = "Параметры :\n{: ^7}{: ^7}{: ^7}{: ^7}{: ^8}\n".format("N1", "N2", "h1", "h2", "eps")
        s += "{: ^7d}{: ^7d}{: ^7.3f}{: ^7.3f}{: ^8}\n\n".format(N1, N2, h1, h2, eps)
        s += "Результаты:\n {: ^12}{: ^24}\n".format("Iter count", "Max Fail")
        s += " {: ^12}{: ^24.15f}\n{}\n".format(iter_count, np.max(yerr), "="*36)
    else:
        s = "{: ^7}{: ^7}{: ^7}{: ^7}{: ^8}".format("N1", "N2", "h1", "h2", "eps")
        s += "{: ^12}{: ^24}\n".format("Iter count", "Max Fail")
        s += "{: ^7d}{: ^7d}{: ^7.3f}{: ^7.3f}{: ^8}".format(N1, N2, h1, h2, eps)
        s += " {: ^12}{: ^24.15f}\n{}\n".format(iter_count, np.max(yerr), "="*72)
    return s

def edge_computing(ys, u, x1, x2):
    N1, N2 = np.shape(ys)
    for i, x in enumerate(x1):
        ys[i, 0] = u(x,0)
        ys[i, -1] = u(x,1)
    for i, x in enumerate(x2):
        ys[0, i] = u(0, x)
        ys[-1, i] = u(1, x)

def final_error(ys, ysol):
    yerr = ys.copy()
    N1, N2 = np.shape(ys)

    for i in range(N1):
        for j in range(N2):
            yerr[i, j] = abs(ysol[i, j]-ys[i, j])
    return yerr

def test_solution(ys, Cs, F):
    rka = 0.0
    N1, N2 = np.shape(ys)
    for i in range(1, N1-1):
        for j in range(1, N2-1):
            rka = max(rka, abs(F[i,j] + Cs[i][j][1]*ys[i+1, j] + \
                      + Cs[i][j][2]*ys[i-1, j] + Cs[i][j][3]*ys[i, j+1] + \
                      + Cs[i][j][4]*ys[i, j-1] - Cs[i][j][0]*ys[i, j]))

    return rka

def my3d_plot(zdata, x1, x2):
    X, Y = np.meshgrid(x1, x2)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, zdata, color='0.75', rstride=1, cstride=1, cmap="magma")
    ax.view_init(elev=10., azim=0)
    plt.show()    
