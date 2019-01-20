#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mymetods.h"
#include "task.h"

#define MPI_PROC_NULL 1

const size_t N1 = 20, N2 = 20;

int main(int argc, char const *argv[])
{
    double *x1 = linspace(0, 1, N1);
    double *x2 = linspace(0, 1, N2);

    const double h1 = 1.0/N1;
    const double h2 = 1.0/N2;
    const double eps = 1e-5;
    const double eps_j = 1e-5;

    const size_t maxiter = 100;
    const size_t maxiter_jacobi = 300;

    int start_row,last_row,num_row,start_col,last_col,num_col;
    int proc_left,proc_right, proc_down,proc_up;
    int limits1[2], limits2[2];

    proc_left = proc_right = proc_down = proc_up = 1;
    /* rows of matrix I have to process */
    start_row = 0;
    last_row = N1-1;
    num_row = last_row - start_row + 1;
    /* columns of matrix I have to process */
    start_col = 0;
    last_col = N2-1;
    num_col = last_col - start_col + 1;

    double ** ys = calloc(num_row+2, sizeof(double*));
    ys[0] = calloc((num_row+2) * (num_col+2), sizeof(double));
    for (int i = 1; i < num_row+2; i++) {
        ys[i] = ys[0] + i * (num_col+2);
    }
    double ** ys1 = calloc(num_row+2, sizeof(double*));
    ys1[0] = calloc((num_row+2) * (num_col+2), sizeof(double));
    for (int i = 1; i < num_row+2; i++) {
        ys1[i] = ys1[0] + i * (num_col+2);
    }

    double** ysol = calloc(num_row+2, sizeof(double*));
    double* Csi;
    double*** Cs = calloc(num_row+2, sizeof(double**));
    double** F = calloc(num_row+2, sizeof(double*));
    for (int i=0; i<num_row+2; i++)
    {
        ysol[i] = calloc(num_col+2, sizeof(double));
        F[i] = calloc(num_col+2, sizeof(double));
        Cs[i] = calloc(num_col+2, sizeof(double*));
        for (int j = 0; j < num_col+2; j++)
            Cs[i][j] = calloc(5, sizeof(double));
    }

    for(int i = 1; i < num_row+1; i++)
        for(int j=1; j < num_col+1; j++)
        {
            if ((i==1)&&(proc_up==MPI_PROC_NULL)) {
                ys[i][j] = u(0, x2[start_col + j - 1]);
                continue;
            }
            if ((i==num_row)&&(proc_down==MPI_PROC_NULL)) {
                ys[i][j] = u(1, x2[start_col + j - 1]);
                continue;
            }
            if ((j==1)&&(proc_left==MPI_PROC_NULL)) {
                ys[i][j] = u(x1[start_row + i - 1], 0);
                continue;
            }
            if ((j==num_col)&&(proc_right==MPI_PROC_NULL)) {
                ys[i][j] = u(x1[start_row + i - 1], 1);
                continue;
            }
        }
    memcpy(ys1[0],ys[0], (num_row+2) * (num_col+2)*sizeof(double));

    int iter_count;
    for (iter_count = 0; iter_count < maxiter;  iter_count++) {
        limits1[0] = (proc_up==MPI_PROC_NULL) ? 2 : 1; limits1[1] = (proc_down==MPI_PROC_NULL) ? num_row : num_row+1;
        limits2[0] = (proc_left==MPI_PROC_NULL) ? 2 : 1; limits2[1] = (proc_right==MPI_PROC_NULL) ? num_col : num_col+1;
        for (int i = limits1[0]; i <  limits1[1]; i++) {
            for (int j = 2; j <  num_col; j++){
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

        if (test_solution(ys, Cs, F, num_row, num_col) < eps) break;


        for (int iter_count_j = 0;  iter_count_j < maxiter_jacobi; iter_count_j++) {
            limits1[0] = (proc_up==MPI_PROC_NULL) ? 2 : 1; limits1[1] = (proc_down==MPI_PROC_NULL) ? num_row : num_row+1;
            limits2[0] = (proc_left==MPI_PROC_NULL) ? 2 : 1; limits2[1] = (proc_right==MPI_PROC_NULL) ? num_col : num_col+1;
            for (int i = limits1[0]; i <  limits1[1]; i++) {
                for (int j = limits2[0]; j <  limits2[1]; j++){
                    ys1[i][j] = (F[i][j] + Cs[i][j][1]*ys[i+1][j] + Cs[i][j][2]*ys[i-1][j] +
                        + Cs[i][j][3]*ys[i][j+1] + Cs[i][j][4]*ys[i][j-1])/Cs[i][j][0];

                }
            }
            memcpy(ys[0],ys1[0], (num_row+2) * (num_col+2)*sizeof(double));
            if (test_solution(ys, Cs, F, num_row, num_col) < eps_j) break;
        }
    }

    for (int i = 0; i <  num_row+2; i++) {
        printf("\n");
        for (int j = 0; j <  num_col+2; j++){
            printf("%6.3f ", ys[i][j]);

        }
    }
    printf("\n");

    //Задаем точное решение
    solution(x1+start_row, num_row, x2+start_col, num_col, ysol);

    //вычисляем ошибку на сетке
    double yerr = final_error(ys, ysol, num_row, num_col);

    // выводим количество итераций и максимальную ошибку
    print_res(num_row, num_col, h1, h2, eps, iter_count, yerr);


    for(int i=0; i<num_row+2; i++)
    {
        for (int j = 0; j < num_col+2; j++)
            free(Cs[i][j]);
        free(Cs[i]);
        free(F[i]);
        free(ysol[i]);
    }
    free(Cs);
    free(ys[0]);
    free(ys);
    free(ysol);
    free(F);
    free(x1);
    free(x2);
    return 0;
}
