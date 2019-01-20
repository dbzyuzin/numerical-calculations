#include <math.h>
#include <stdio.h>
#include "mymetods.h"
#include "task.h"
#include "mpi.h"

void grid(int nx, int ny, int np, int* o_np1, int* o_np2)
{
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

int MyNetInit(int* argc, char*** argv, int* np, int* mp)
{
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
    for (int i=0; i<N; i++){
        x[i] = xa + i*h;
    }
	return x;
}

double dmax(const double x1, const double x2)
{
	if (x1 > x2)
		return x1;
	return x2;
}

double test_solution(double** ys, double*** Cs, double** F, const size_t n1, const size_t n2)
{
    double rka = 0.0;

    for (int i = 1; i < n1+1; i++){
        for (int j = 1; j < n2+1; j++){
            rka = dmax(rka, fabs(F[i][j] + Cs[i][j][1]*ys[i+1][j] +
                      + Cs[i][j][2]*ys[i-1][j] + Cs[i][j][3]*ys[i][j+1] +
                      + Cs[i][j][4]*ys[i][j-1] - Cs[i][j][0]*ys[i][j]));
		}
	}

    return rka;
}

void print_res(const int n1, const int n2,
                const double h1, const double h2,
                const double eps, const int iter_count, const double rka)
{
    printf("Параметры :\n%7s %7s %7s %7s %8s\n","N1", "N2", "h1", "h2", "eps");
	printf("%7d %7d %7.3f %7.3f %8.0e\n\n",n1, n2, h1, h2, eps);
	printf("Результаты:\n %10s %24s \n", "Iter count", "Max Fail");
	printf(" %10d %24.10f\n", iter_count, rka);
}
void fprint_res(const int n1, const int n2,
                const double h1, const double h2,
                const double eps, const int iter_count, const double rka)
{
	FILE *f = fopen("test.txt", "a");
	if (f == NULL){
	    printf("Error opening file!\n");
		return;
	}
	fprintf(f, "Параметры :\n%7s %7s %7s %7s %8s\n","N1", "N2", "h1", "h2", "eps");
	fprintf(f, "%7d %7d %7.3f %7.3f %8.0e\n\n",n1, n2, h1, h2, eps);
	fprintf(f, "Результаты:\n %10s %24s \n", "Iter count", "Max Fail");
	fprintf(f, " %10d %24.10f\n", iter_count, rka);

}

void solution(double* x1, const size_t n1,
            double* x2, const size_t n2, double** ysol)
{
    for (int i=1; i<n1+1; i++)
        for (int j=1; j<n2+1; j++)
            ysol[i][j] = u(x1[i-1],x2[j-1]);
}

//not for parallel version
void edge_computing(double* x1, const size_t n1,
             double* x2, const size_t n2, double** ys)
{
    for(int i=1; i < n1+1; i++) {
        ys[i][1] = u(x1[i-1],0);
        ys[i][n2] = u(x1[i-1],1);
    }
    for(int i=1; i < n2+1; i++) {
        ys[1][i] = u(0, x2[i-1]);
        ys[n1][i] = u(1, x2[i-1]);
    }

}

double final_error(double** ys, double** ysol, const size_t n1, const size_t n2)
{
    double err = 0;
    for (int i=1; i<n1+1; i++)
        for (int j=1; j<n2+1; j++)
            err = dmax(fabs(ysol[i][j]-ys[i][j]), err);
    return err;
}
