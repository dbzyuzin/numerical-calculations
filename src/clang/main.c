#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mymetods.h"
#include "task.h"
#include "mpi.h"
#include "omp.h"

int main(int argc, char* argv[])
{
    int mp, np;

    const size_t N1 = 200, N2 = 200;

    double* x1 = linspace(0, 1, N1);
    double* x2 = linspace(0, 1, N2);

    const double h1 = 1.0/N1;
    const double h2 = 1.0/N2;
    const double eps = 1e-5;

    const size_t maxiter = 100;
    const size_t maxiter_jacobi = 300;

    int start_row,last_row,num_row,start_col,last_col,num_col;
    int proc_left,proc_right, proc_down,proc_up;
    int limits1[2], limits2[2];
    double t1;

    int isper[] = {0,0}; //периодичность решетки
    int dim[] = {0,0};         //размерность
    int coords[2];
    MPI_Request req[8];
    MPI_Status status[8];
    MPI_Comm newcomm;
    MPI_Datatype vectype;

    int ret = MyNetInit(&argc, &argv, &np, &mp);
    if (ret)
        return -1;

    grid(N1, N2, np, &dim[0], &dim[1]);
    MPI_Cart_create(MPI_COMM_WORLD,2,dim,isper,1,&newcomm);
    MPI_Cart_shift(newcomm,0,1,&proc_up,&proc_down);
    MPI_Cart_shift(newcomm,1,1,&proc_left,&proc_right);
    MPI_Comm_rank (newcomm, &mp); /* my place in MPI system */
    MPI_Cart_coords(newcomm,mp,2,coords);

    /* rows of matrix I have to process */
    start_row = (coords[0] * N1) / dim[0];
    last_row = (((coords[0] + 1) * N1) / dim[0])-1;
    num_row = last_row - start_row + 1;
    /* columns of matrix I have to process */
    start_col = (coords[1] * N2) / dim[1];
    last_col = (((coords[1] + 1) * N2) / dim[1])-1;
    num_col = last_col - start_col + 1;

    MPI_Type_vector(num_row,1,num_col+2,MPI_DOUBLE,&vectype);
    MPI_Type_commit(&vectype);

    double** ys = calloc(num_row+2, sizeof(double*));
    ys[0] = calloc((num_row+2) * (num_col+2), sizeof(double));
    for (size_t i = 1; i < num_row+2; i++)
        ys[i] = ys[0] + i * (num_col+2);

    double** ys1 = calloc(num_row+2, sizeof(double*));
    ys1[0] = calloc((num_row+2) * (num_col+2), sizeof(double));
    for (size_t i = 1; i < num_row+2; i++)
        ys1[i] = ys1[0] + i * (num_col+2);

    double** ysol = calloc(num_row+2, sizeof(double*));
    double* Csi;
    double*** Cs = calloc(num_row+2, sizeof(double**));
    double** F = calloc(num_row+2, sizeof(double*));
    for (size_t i=0; i<num_row+2; i++)
    {
        ysol[i] = calloc(num_col+2, sizeof(double));
        F[i] = calloc(num_col+2, sizeof(double));
        Cs[i] = calloc(num_col+2, sizeof(double*));
        for (int j = 0; j < num_col+2; j++)
            Cs[i][j] = calloc(5, sizeof(double));
    }

    if (mp==0)
        printf("CALC STARTED on %d*%d processors with %d*%d array, my part %dx%d\n",
            dim[0], dim[1], N1, N2, num_row, num_col);

    for(size_t i = 1; i < num_row+1; i++)
        for(size_t j=1; j < num_col+1; j++)
        {
            if ((i==1)&&(proc_up==MPI_PROC_NULL)) 
            {
                ys[i][j] = u(0, x2[start_col + j - 1]);
                continue;
            }
            if ((i==num_row)&&(proc_down==MPI_PROC_NULL)) 
            {
                ys[i][j] = u(1, x2[start_col + j - 1]);
                continue;
            }
            if ((j==1)&&(proc_left==MPI_PROC_NULL)) 
            {
                ys[i][j] = u(x1[start_row + i - 1], 0);
                continue;
            }
            if ((j==num_col)&&(proc_right==MPI_PROC_NULL)) 
            {
                ys[i][j] = u(x1[start_row + i - 1], 1);
                continue;
            }
        }
    memcpy(ys1[0], ys[0], (num_row+2)*(num_col+2)*sizeof(double));

    limits1[0] = (proc_up==MPI_PROC_NULL) ? 2 : 1; 
    limits1[1] = (proc_down==MPI_PROC_NULL) ? num_row : num_row+1;
    limits2[0] = (proc_left==MPI_PROC_NULL) ? 2 : 1;
    limits2[1] = (proc_right==MPI_PROC_NULL) ? num_col : num_col+1;

    // iteration loop
    MPI_Barrier(newcomm);
    t1 = MPI_Wtime();
    size_t iter_count;
    for (iter_count = 0; iter_count < maxiter;  iter_count++) 
    {
        for (size_t i = limits1[0]; i < limits1[1]; i++)
            for (size_t j = limits2[0]; j < limits2[1]; j++)
            {
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

        if (test_solution(ys, Cs, F, num_row, num_col) < eps) 
            break;

        for (size_t iter_count_j = 0;  iter_count_j < maxiter_jacobi; iter_count_j++) 
        {
            for (size_t i = limits1[0]; i <  limits1[1]; i++) {
                for (size_t j = limits2[0]; j <  limits2[1]; j++){
                    ys1[i][j] = (F[i][j] + Cs[i][j][1]*ys[i+1][j] + Cs[i][j][2]*ys[i-1][j] +
                        + Cs[i][j][3]*ys[i][j+1] + Cs[i][j][4]*ys[i][j-1])/Cs[i][j][0];
                }
            }

            memcpy(ys[0],ys1[0], (num_row+2) * (num_col+2)*sizeof(double));
            MPI_Irecv(&ys[0][1],num_col,MPI_DOUBLE,
            proc_up, 1215, MPI_COMM_WORLD, &req[0]);
            MPI_Isend(&ys[num_row][1],num_col,MPI_DOUBLE,
            proc_down, 1215, MPI_COMM_WORLD,&req[1]);
            MPI_Irecv(&ys[num_row+1][1],num_col,MPI_DOUBLE,
            proc_down, 1216, MPI_COMM_WORLD, &req[2]);
            MPI_Isend(&ys[1][1],num_col,MPI_DOUBLE,
            proc_up, 1216, MPI_COMM_WORLD,&req[3]);
            MPI_Irecv(&ys[1][0],1,vectype,
            proc_left, 1217, MPI_COMM_WORLD, &req[4]);
            MPI_Isend(&ys[1][num_col],1,vectype,
            proc_right, 1217, MPI_COMM_WORLD,&req[5]);
            MPI_Irecv(&ys[1][num_col+1],1,vectype,
            proc_right, 1218, MPI_COMM_WORLD, &req[6]);
            MPI_Isend(&ys[1][1],1,vectype,
            proc_left, 1218, MPI_COMM_WORLD,&req[7]);
            MPI_Waitall(8,req,status);
        }
    }

    //Задаем точное решение
    solution(x1+start_row, num_row, x2+start_col, num_col, ysol);

    //вычисляем ошибку на сетке
    double yerr = final_error(ys, ysol, num_row, num_col);
    if (np>1)
        double yerr1 = yerr; MPI_Allreduce(&yerr1,&yerr,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    //Вычисляем общее время работы
    double time = MPI_Wtime()-t1;
    if (np>1)
        double time1 = time; MPI_Allreduce(&time1,&time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    // выводим количество итераций и максимальную ошибку
    if (mp==0) 
        print_res(N1, N2, h1, h2, eps, iter_count, yerr, time);

    for(size_t i=0; i<num_row+2; i++)
    {
        for (size_t j = 0; j < num_col+2; j++)
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
    MPI_Finalize();
    return 0;
}
