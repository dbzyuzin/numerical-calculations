#include "task.h"

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
    return M_05_PI2*(1-u*u);
}

double u(const double x1, const double x2)
{
    return cos(M_05_PI*(x1+x2));
}

double f(const double u)
{
    return M_PI2*u*u*u;
}
