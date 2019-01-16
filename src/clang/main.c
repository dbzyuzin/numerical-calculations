#include <stdbool.h>
#include <stdlib.h>
#include "mymetods.h"
#include "task.h"

int main(int argc, char const *argv[])
{
    const size_t N1 = 20, N2 = 20;
    const double *restrict x1 = linspace(0, 1, N1);
    const double *restrict x2 = linspace(0, 1, N2);

    const double h1 = 1.0/N1;
    const double h2 = 1.0/N2;
    const double eps = 1e-5;
    const double eps_j = 1e-5;


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
        Cs[i] = calloc(N2, sizeof(double*));
        for (int j = 0; j < N2; j++)
            Cs[i][j] = calloc(5, sizeof(double));
    }
    edge_computing(x1, N1, x2, N2, ys);

    int iter_count;
    for (iter_count = 0; iter_count < maxiter;  iter_count++) {

        for (int i = 1; i <  N1-1; i++) {
            for (int j = 1; i <  N2-1; i++){
                Csi = Cs[i][j];
                Csi[0] = (h2/h1*(ki(ys[i+1][j], ys[i][j]) + ki(ys[i-1][j], ys[i][j])) +
                          + h1/h2*(ki(ys[i][j+1], ys[i][j]) +
                          + ki(ys[i][j-1], ys[i][j])) +
                          + h1*h2*q(ys[i][j]));
                Csi[1] = h2/h1*ki(ys[i+1][j], ys[i][j]);
                Csi[2] = h2/h1*ki(ys[i-1][j], ys[i][j]);
                Csi[3] = h1/h2*ki(ys[i][j+1], ys[i][j]);
                Csi[4] = h1/h2*ki(ys[i][j-1], ys[i][j]);

                F[i][j] = h1*h2*f(ys[i][j]);
                Cs[i][j] = Csi;
            }
        }

        if (test_solution(ys, Cs, F, N1, N2) < eps) break;


        for (int iter_count_j = 0;  iter_count_j < maxiter_jacobi; iter_count_j++) {
            for (int i = 1; i <  N1-1; i++) {
                for (int j = 1; i <  N2-1; i++){
                    ys[i][j] = (F[i][j] + Cs[i][j][1]*ys[i+1][j] + Cs[i][j][2]*ys[i-1][j] +
                        + Cs[i][j][3]*ys[i][j+1] + Cs[i][j][4]*ys[i][j-1])/Cs[i][j][0];

            if (test_solution(ys, Cs, F, N1, N2) < eps_j) break;
                }
            }
        }
    }


    //Задаем точное решение
    solution(x1, N1, x2, N2, ysol);

    //вычисляем ошибку на сетке
    double yerr = final_error(ys, ysol, N1, N2);

    // выводим количество итераций и максимальную ошибку
    print_res(N1, N2, h1, h2, eps, iter_count, yerr);


    for(int i=0; i<N1; i++)
    {
        for (int j = 0; j < N2; j++)
            free(Cs[i][j]);
        free(Cs[i]);
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
