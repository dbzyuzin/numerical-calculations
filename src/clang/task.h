#ifndef TASK
#define TASK
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define M_PI2 (M_PI*M_PI)
#define M_05_PI2 (0.5*M_PI2)
#define M_05_PI (0.5*M_PI)

double k(const double u);
double ki(const double u1, const double u2);

double q(const double u);

double u(const double x1, const double x2);

double f(const double u);
#endif
