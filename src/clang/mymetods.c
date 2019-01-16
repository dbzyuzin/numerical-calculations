
#include <stdlib.h>
#include "mymetods.h"

double* linspace(int xa, int xb, int N) {
	double *x = calloc(N, sizeof(double));
    double h = (xb - xa)/(double)(N-1);
    for (int i=0; i<N; i++){
        x[i] = xa + i*h;
    }
	return x;
}



void print_res(int N1, int N2, double h1, double h2, double eps, int iter_count, double rka) {
    printf("Параметры :\n\n%7s %7s %7s %7s %8s\n","N1", "N2", "h1", "h2", "eps");
	printf("%7d %7d %7.3f %7.3f %8.0e\n\n",N1, N2, h1, h2, eps);
	printf("Результаты:\n\n %10s %24s \n", "Iter count", "Max Fail");
	printf(" %10d %24.10f\n", iter_count, rka);
}

void solution(const double* restrict x1, const size_t N1,
            const double* restrict x2, const size_t N2, double* ysol)
{
    for (int i=0; i<N1; i++)
        for (int j=0; j<N2; j++)
            ysol[i*N1+j] = u(x1[i],x2[j]);
}

void edge_computing(const double* restrict x1, const size_t N1,
            const double* restrict x2, const size_t N2)
{
    for(int i=0; i < N1; i++) {
        ys[i*N1] = u(x1[i],0);
        ys[i*N1+N2-1] = u(x1[i],1);
    }
    for(int i=0; i < N2; i++) {
        ys[i*N2] = u(0, x2[i]);
        ys[N1-1+i*N2] = u(1, x2[i]);
    }
}
