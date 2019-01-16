#include <math.h>
#define M_PI 3.14159265358979323846
#define M_PI2 (M_PI*M_PI)
#define M_05_PI2 (0.5*M_PI2)
#define M_05_PI (0.5*M_PI)

inline double k(const double u) {
    return 1.0 + u*u;
}
double ki(const double u1, const double u2) {
    const double k1 = k(u1);
    const double k2 = k(u2);
    return 2*(k1*k2)/(k1+k2);
}

inline double q(const double u) {
    return M_05_PI2*(1-u*u);
}
double u(const double x1, const double x2) {
    return cos(M_05_PI*(x1+x2));
}

inline double f(const double u) {
    return M_PI2*u*u*u;
}

void solution(double (*u)(const double, const double), double* x1, size_t N1, 
            double* x2, size_t N2, double* ysol) 
{
    for (int i=0; i<N1; i++)
        for (int j=0; j<N2; j++)
            ysol[i*N1+j] = u(x1[i],x2[j]);
}
