#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{
    const size_t N1 = 20, N2 = 20;
    const double *restrict x1 = linspace(0, 1, N1);
    const double *restrict x2 = linspace(0, 1, N2);

    const double h1 = 1.0/N1;
    const double h2 = 1.0/N2;
    const double eps = 1e-5;

    const size_t maxiter = 100;
    const size_t maxiter_jacobi = 333;

    double** ys = calloc(N1, sizeof(double*));
    double** ysol = calloc(N1, sizeof(double*));
    double* Csi;
    double*** Cs = calloc(N1, sizeof(double**));
    double** F = calloc(N1, sizeof(double*));
    for (int i=0; i<N1; i++) 
    {
        ys[i] = calloc(N2, sizeof(double));
        ysol[i] = calloc(N2, sizeof(double));
        F[i] = calloc(N2, sizeof(double));
        ys[i] = calloc(N2, sizeof(double*));
        for (int j = 0; j < N2; j++)
            ys[i][j] = calloc(5, sizeof(double));
    }
    edge_computing(x1, N1, x2, N2, ys);






    //Задаем точное решение
    solution(x1, N1, x2, N2, ysol);

    //вычисляем ошибку на сетке
    yerr = final_error(ys, ysol, N1, N2);

    // выводим количество итераций и максимальную ошибку
    print_res(N1, N2, h1, h2, eps, iter_count, yerr)


    for(int i=0; i<N1; i++) 
    {
        for (int j = 0; j < N2; j++)
            free(ys[i][j]);
        free(ys[i]);
        free(F[i]);
        free(ys[i]);
        free(ysol[i]);
    }
    free(Cs);
    free(ys);
    free(ysol);
    free(F);
    free(x1);
    free(x2);
    return 0;
}
