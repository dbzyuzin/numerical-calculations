#ifndef MYMETHODS
#define MYMETHODS
#include <stdlib.h>

double* linspace(int xa, int xb, int N);

double dmax(const double x1, const double x2);

double test_solution(double** ys, double*** Cs, double** F, const size_t n1, const size_t n2);
void print_res(const int n1, const int n2,
                const double h1, const double h2,
                const double eps, const int iter_count, const double rka);
void fprint_res(const int n1, const int n2,
                const double h1, const double h2,
                const double eps, const int iter_count, const double rka);

void solution(double* x1, const size_t n1,
            double* x2, const size_t n2, double** ysol);

void edge_computing(double* x1, const size_t n1,
             double* x2, const size_t n2, double** ys);

double final_error(double** ys, double** ysol, const size_t n1, const size_t n2);

int MyNetInit(int* argc, char*** argv, int* np, int* mp);

void grid(int nx, int ny, int np, int* o_np1, int* o_np2);
#endif
