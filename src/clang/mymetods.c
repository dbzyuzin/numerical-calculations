#include <stdlib.h>
#include <math.h>
#include "mymetods.h"

double* linspace(int xa, int xb, int N) {
	double *x = calloc(N, sizeof(double));
    double h = (xb - xa)/(double)(N-1);
    for (int i=0; i<N; i++){
        x[i] = xa + i*h;
    }
	return x;
}

double test_solution(double **ys, double ***Cs, double **F,int N1,int N2) {
    double rka = 0.0

    for (int i = 1; i < N1-1; i++){
        for (int j = 1; i < N2-1; i++){
            rka = dmax(rka, abs(F[i][j] + Cs[i][j][1]*ys[i+1][j] + \
                      + Cs[i][j][2]*ys[i-1][j] + Cs[i][j][3]*ys[i][j+1] + \
                      + Cs[i][j][4]*ys[i][j-1] - Cs[i][j][0]*ys[i][j]))
		}
	}

    return rka;
}

double dmax(const double x1, const double x2) {
	if (x1 > x2)
		return x1;
	return x2;
}

void print_res(const int N1, const int N2, 
                const double h1, const double h2,
                const double eps, const int iter_count, const double rka) 
{
    printf("Параметры :\n\n%7s %7s %7s %7s %8s\n","N1", "N2", "h1", "h2", "eps");
	printf("%7d %7d %7.3f %7.3f %8.0e\n\n",N1, N2, h1, h2, eps);
	printf("Результаты:\n\n %10s %24s \n", "Iter count", "Max Fail");
	printf(" %10d %24.10f\n", iter_count, rka);
}

void solution(const double* restrict x1, const size_t N1,
            const double* restrict x2, const size_t N2, double** ysol)
{
    for (int i=0; i<N1; i++)
        for (int j=0; j<N2; j++)
            ysol[i][j] = u(x1[i],x2[j]);
}

void edge_computing(const double* restrict x1, const size_t N1,
            const double* restrict x2, const size_t N2, double**restrict ys)
{
    for(int i=0; i < N1; i++) {
        ys[i][0] = u(x1[i],0);
        ys[i][N2-1] = u(x1[i],1);
    }
    for(int i=0; i < N2; i++) {
        ys[0][i] = u(0, x2[i]);
        ys[N1-1][i] = u(1, x2[i]);
    }
}

double final_error(const double* ys, const double* ysol, const size_t N1, const size_t N2)
{
    double err = 0;
    for (int i=0; i<N1; i++)
        for (int j=0; j<N2; j++)
            err = dmax(abs(ysol[i, j]-ys[i, j]), err);
    return err;
}