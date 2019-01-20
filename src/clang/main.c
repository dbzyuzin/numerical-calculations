/*Jacobi-2d program */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#define m_printf if (mp==0)printf
#define M_PI 3.14159265358979323846
#define M_05_PI (0.5*M_PI)
int i,j,it;

void grid(int nx, int ny, int np, int* o_np1, int* o_np2) {
    double s;
    int np1, np2;
    if (np==1) {
        np1 = 1;
        np2 = 1;
    } else {
      s = sqrt(((double)np)) * ((double)nx) / ((double)ny);
      np1 = floor(s);
      if (s>0.5+((double)np1)) np1++;
      np2 = np / np1;
      if (np1*np2!=np) {
        if (nx>ny) {np1 = np; np2 = 1;} else {np1 = 1; np2 = np;}
      }
    }
    *o_np1 = np1;
    *o_np2 = np2;
}

int MyNetInit(int* argc, char*** argv, int* np, int* mp) {
    int i;

    i = MPI_Init(argc,argv);
    if (i != 0){
        fprintf(stderr,"MPI initialization error");
        MPI_Finalize();
        exit(i);
    }

    MPI_Comm_size(MPI_COMM_WORLD,np);
    MPI_Comm_rank(MPI_COMM_WORLD,mp);

    // sleep(1);

    return 0;
}

double dmax(const double x1, const double x2)
{
	if (x1 > x2)
		return x1;
	return x2;
}

double* linspace(int xa, int xb, int N)
{
	double *x = calloc(N, sizeof(double));
    double h = (xb - xa)/(double)(N-1);
    int i;
    for (i=0; i<N; i++){
        x[i] = xa + i*h;
    }
	return x;
}
double u(const double x1, const double x2)
{
    return cos(M_05_PI*(x1+x2));
}

double k(const double u)
{
    return 1.0 + u*u;
}

double ki(const double u1, const double u2)
{
    const double k1 = k(u1);
    const double k2 = k(u2);
    return 2*(k1*k2)/(k1+k2);
}

double q(const double u)
{
    return M_05_PI*M_05_PI*(1-u*u);
}

double f(const double u)
{
    return M_PI*M_PI*u*u*u;
}

void solution(double*  x1, const size_t num_row,
            double* x2, const size_t num_col, double** ysol)
{
    for(i=1; i<=num_row; i++)
        for(j=1; j<=num_col; j++)
            ysol[i][j] = u(x1[i-1],x2[j-1]);
}



int main(int argc, char **argv)
{
    MPI_Request req[8];
    MPI_Status status[8];
    MPI_Comm newcomm;
    MPI_Datatype vectype;

    int mp, np;
    int start_row,last_row,num_row,start_col,last_col,num_col;
    int proc_left,proc_right, proc_down,proc_up;
    int isper[] = {0,0}; //периодичность решетки
    int dim[2]; //размерность
    int coords[2];

    double t1;
    int np1, np2;

    const size_t N1 = 100, N2 = 100;

    const double h1 = 1.0/N1;
    const double h2 = 1.0/N2;
    const double eps = 1e-5;
    const double eps_j = 1e-5;


    const size_t maxiter = 500;
    const size_t maxiter_jacobi = 300;

    MyNetInit(&argc, &argv, &np, &mp);
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

    double (* A)[num_col+2];
    double (* B)[num_col];
    double *x1 = linspace(0, 1, N1);
    double *x2 = linspace(0, 1, N2);

    double (* ys)[num_col+2];
    ys = malloc ((num_row+2) * (num_col+2) * sizeof(double));

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

    m_printf("JAC2 STARTED on %d*%d processors with %d*%d array, it=%d\n",dim[0],dim[1],N1,N2,maxiter);


    /* dynamically allocate data structures */
    A = malloc ((num_row+2) * (num_col+2) * sizeof(double));
    B = malloc (num_row * num_col * sizeof(double));

    for(i=1; i<=num_row; i++)
        for(j=1; j<=num_col; j++)
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

    /****** iteration loop *************************/
    MPI_Barrier(newcomm);
    t1=MPI_Wtime();
    for(it=1; it<=maxiter; it++)
    {
        for(i=1; i<=num_row; i++)
        {
            if (((i==1)&&(proc_up==MPI_PROC_NULL))||((i==num_row)&&(proc_down==MPI_PROC_NULL))) continue;
            for(j=1; j<=num_col; j++)
            {
                if (((j==1)&&(proc_left==MPI_PROC_NULL))||((j==num_col)&&(proc_right==MPI_PROC_NULL))) continue;
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

        for (int it1 = 0;  it1 < maxiter_jacobi; it1++) {
            for(i=1; i<=num_row; i++)
            {
                if (((i==1)&&(proc_up==MPI_PROC_NULL))||((i==num_row)&&(proc_down==MPI_PROC_NULL))) continue;
                for(j=1; j<=num_col; j++)
                {
                    if (((j==1)&&(proc_left==MPI_PROC_NULL))||((j==num_col)&&(proc_right==MPI_PROC_NULL))) continue;
                    ys[i][j] = (F[i][j] + Cs[i][j][1]*ys[i+1][j] + Cs[i][j][2]*ys[i-1][j] +
                        + Cs[i][j][3]*ys[i][j+1] + Cs[i][j][4]*ys[i][j-1])/Cs[i][j][0];
                }
            }
        }

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

        // if (coords[0]==1 && coords[1]==1 && it==1)
        //     for(i=0; i<=num_row+1; i++)
        //     {
        //         printf("\n");
        //         for(j=0; j<=num_col+1; j++)
        //         {
        //             printf("%5.3f ", ys[i][j]);
        //         }
        //     }
        // printf("\n");

    }

    //Задаем точное решение
    solution(x1+start_row, num_row, x2+start_col, num_col, ysol);

    //вычисляем ошибку на сетке
    double err = 0;
    for(i=1; i<=num_row; i++)
        for(j=1; j<=num_col; j++)
            // if (fabs(ysol[i][j]-ys[i][j]) > err) {
            //     err = fabs(ysol[i][j]-ys[i][j]);
            //     if (coords[0]==1 && coords[1]==1)
            //         printf("err: %f,      %dx%d,      %fx%f \n", err, i, j, ysol[i][j], ys[i][j]);
            // }
            err = dmax(fabs(ysol[i][j]-ys[i][j]), err);

    printf("%d: Time of task=%lf, err %f, %d\n",mp,MPI_Wtime()-t1, err, it);
    for(int i=0; i<num_row+2; i++)
    {
        for (int j = 0; j < num_col+2; j++)
            free(Cs[i][j]);
        free(Cs[i]);
        free(F[i]);
        free(ysol[i]);
    }
    free(Cs);
    free(ys);
    free(ysol);
    free(F);
    free(x1);
    free(x2);
    MPI_Finalize();
    return 0;
}
