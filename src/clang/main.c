/*Jacobi-2d program */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
#define m_printf if (mp==0)printf
#define L 6
#define LC 2
#define ITMAX 100
int i,j,it,k;

int main(int argc, char **argv)
{
    double (* A)[L/LC+2];
    double (* B)[L/LC];

    MPI_Request req[8];
    int mp, np;
    int start_row,last_row,num_row,start_col,last_col,num_col;
    MPI_Status status[8];
    double t1;
    int isper[] = {0,0}; //периодичность решетки
    int dim[2]; //размерность
    int coords[2];
    MPI_Comm newcomm;
    MPI_Datatype vectype;
    int proc_left,proc_right, proc_down,proc_up;
    int np1, np2, nx, ny;
    double s;

    nx = ny = L;

    MPI_Init (&argc, &argv); /* initialize MPI system */
    MPI_Comm_size (MPI_COMM_WORLD, &np); /* size of MPI system */
    MPI_Comm_rank (MPI_COMM_WORLD, &mp); /* my place in MPI system */

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

    dim[0]=np/LC;
    dim[1]=LC;
    if ((L%dim[0])||(L%dim[1]))
    {
        m_printf("ERROR: array[%d*%d] is not distributed on %d*%d processors\n",L,L,dim[0],dim[1]);
        MPI_Finalize();
        exit(1);
    }
    MPI_Cart_create(MPI_COMM_WORLD,2,dim,isper,1,&newcomm);
    MPI_Cart_shift(newcomm,0,1,&proc_up,&proc_down);
    MPI_Cart_shift(newcomm,1,1,&proc_left,&proc_right);
    MPI_Comm_rank (newcomm, &mp); /* my place in MPI system */
    MPI_Cart_coords(newcomm,mp,2,coords);

    /* rows of matrix I have to process */
    start_row = (coords[0] * L) / dim[0];
    last_row = (((coords[0] + 1) * L) / dim[0])-1;
    num_row = last_row - start_row + 1;
    /* columns of matrix I have to process */
    start_col = (coords[1] * L) / dim[1];
    last_col = (((coords[1] + 1) * L) / dim[1])-1;
    num_col = last_col - start_col + 1;
    MPI_Type_vector(num_row,1,num_col+2,MPI_DOUBLE,&vectype);
    MPI_Type_commit(&vectype);
    m_printf("JAC2 STARTED on %d*%d processors with %d*%d array, it=%d\n",dim[0],dim[1],L,L,ITMAX);
    /* dynamically allocate data structures */
    A = malloc ((num_row+2) * (num_col+2) * sizeof(double));
    B = malloc (num_row * num_col * sizeof(double));

    for(i=0; i<=num_row-1; i++)
    {
        for(j=0; j<=num_col-1; j++)
        {
            A[i+1][j+1]=0.;
            B[i][j]=1.+start_row+i+start_col+j;
        }
    }

    /****** iteration loop *************************/
    MPI_Barrier(newcomm);
    t1=MPI_Wtime();
    for(it=1; it<=ITMAX; it++)
    {
        for(i=0; i<=num_row-1; i++)
        {
            if (((i==0)&&(proc_up==MPI_PROC_NULL))||((i==num_row-1)&&(proc_down==MPI_PROC_NULL))) continue;
            for(j=0; j<=num_col-1; j++)
            {
                if (((j==0)&&(proc_left==MPI_PROC_NULL))||((j==num_col-1)&&(proc_right==MPI_PROC_NULL)))
                continue;
                A[i+1][j+1] = B[i][j];
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
                B[i-1][j-1] = (A[i-1][j]+A[i+1][j]+A[i][j-1]+A[i][j+1])/4.;
            }
        }
    }
    printf("%d: Time of task=%lf\n",mp,MPI_Wtime()-t1);
    MPI_Finalize ();
    return 0;
}
