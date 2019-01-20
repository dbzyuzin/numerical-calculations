/*Jacobi-2d program */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#define m_printf if (mp==0)printf
#define M_PI 3.14159265358979323846
#define M_05_PI (0.5*M_PI)
int i,j,it,k;

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
// void edge_computing(const double* restrict x1, const size_t N1,
//             const double* restrict x2, const size_t N2, double** restrict ys)
// {
//     for(int i=0; i < N1; i++) {
//         ys[i][0] = u(x1[i],0);
//         ys[i][N2-1] = u(x1[i],1);
//     }
//     for(int i=0; i < N2; i++) {
//         ys[0][i] = u(0, x2[i]);
//         ys[N1-1][i] = u(1, x2[i]);
//     }
// 	for (int i = 1; i <  N1-1; i++)
//         for (int j = 1; j <  N2-1; j++)
//             ys[i][j] = 0;
//
// }

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

    const size_t N1 = 20, N2 = 20;

    const double h1 = 1.0/N1;
    const double h2 = 1.0/N2;
    const double eps = 1e-5;
    const double eps_j = 1e-5;


    const size_t maxiter = 100;
    const size_t maxiter_jacobi = 20;

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

    double (* A)[num_col+2];
    double (* B)[num_col];
    const double *restrict x1 = linspace(0, 1, N1);
    const double *restrict x2 = linspace(0, 1, N2);

    double** ys = calloc(num_row+2, sizeof(double*));
    for (int i=0; i<num_row+2; i++)
    {
        ys[i] = calloc(num_col+2, sizeof(double));
    }

    MPI_Type_vector(num_row,1,num_col+2,MPI_DOUBLE,&vectype);
    MPI_Type_commit(&vectype);
    m_printf("JAC2 STARTED on %d*%d processors with %d*%d array, it=%d\n",dim[0],dim[1],N1,N2,maxiter);


    /* dynamically allocate data structures */
    A = malloc ((num_row+2) * (num_col+2) * sizeof(double));
    B = malloc (num_row * num_col * sizeof(double));

    for(i=1; i<=num_row; i++)
    {

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
            ys[i][j] = 0;
            // A[i+1][j+1]=0.;
            // B[i][j]=1.+start_row+i+start_col+j;
        }
    }
    if (mp==0)
        for(i=0; i<=num_row+1; i++)
        {
            printf("\n");
            for(j=0; j<=num_col+1; j++)
            {
                printf("%5.3f ", ys[i][j]);
            }
        }

    // edge_computing(x1, num_row, x2, num_col, ys);

    /****** iteration loop *************************/
    MPI_Barrier(newcomm);
    t1=MPI_Wtime();
    for(it=1; it<=maxiter; it++)
    {
        for(i=0; i<=num_row-1; i++)
        {
            if (((i==0)&&(proc_up==MPI_PROC_NULL))||((i==num_row-1)&&(proc_down==MPI_PROC_NULL))) continue;
            for(j=0; j<=num_col-1; j++)
            {
                if (((j==0)&&(proc_left==MPI_PROC_NULL))||((j==num_col-1)&&(proc_right==MPI_PROC_NULL)))
                    continue;
                // A[i+1][j+1] = B[i][j];
            }
        }
        MPI_Irecv(&A[0][1],num_col,MPI_DOUBLE,
        proc_up, 1215, MPI_COMM_WORLD, &req[0]);
        MPI_Isend(&A[num_row][1],num_col,MPI_DOUBLE,
        proc_down, 1215, MPI_COMM_WORLD,&req[1]);
        MPI_Irecv(&A[num_row+1][1],num_col,MPI_DOUBLE,
        proc_down, 1216, MPI_COMM_WORLD, &req[2]);
        MPI_Isend(&A[1][1],num_col,MPI_DOUBLE,
        proc_up, 1216, MPI_COMM_WORLD,&req[3]);
        MPI_Irecv(&A[1][0],1,vectype,
        proc_left, 1217, MPI_COMM_WORLD, &req[4]);
        MPI_Isend(&A[1][num_col],1,vectype,
        proc_right, 1217, MPI_COMM_WORLD,&req[5]);
        MPI_Irecv(&A[1][num_col+1],1,vectype,
        proc_right, 1218, MPI_COMM_WORLD, &req[6]);
        MPI_Isend(&A[1][1],1,vectype,
        proc_left, 1218, MPI_COMM_WORLD,&req[7]);
        MPI_Waitall(8,req,status);

        for(i=1; i<=num_row; i++)
        {
            if (((i==1)&&(proc_up==MPI_PROC_NULL))||
            ((i==num_row)&&(proc_down==MPI_PROC_NULL))) continue;
            for(j=1; j<=num_col; j++)
            {
                if (((j==1)&&(proc_left==MPI_PROC_NULL))||
                ((j==num_col)&&(proc_right==MPI_PROC_NULL))) continue;
                // B[i-1][j-1] = (A[i-1][j]+A[i+1][j]+A[i][j-1]+A[i][j+1])/4.;
            }
        }
    }
    printf("%d: Time of task=%lf\n",mp,MPI_Wtime()-t1);
    MPI_Finalize ();
    return 0;
}
