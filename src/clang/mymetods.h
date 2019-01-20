#ifndef MYMETHODS
#define MYMETHODS
#include <stdlib.h>

double* linspace(int xa, int xb, int N);

double dmax(const double x1, const double x2);

double test_solution(double** ys, double*** Cs, double** F, const size_t N1, const size_t N2);
void print_res(const int N1, const int N2,
                const double h1, const double h2,
                const double eps, const int iter_count, const double rka);
void fprint_res(const int N1, const int N2,
                const double h1, const double h2,
                const double eps, const int iter_count, const double rka);

void solution(double* x1, const size_t N1,
            double* x2, const size_t N2, double** ysol);

void edge_computing(double* x1, const size_t N1,
             double* x2, const size_t N2, double** ys);

double final_error(double** ys, double** ysol, const size_t N1, const size_t N2);
#endif
