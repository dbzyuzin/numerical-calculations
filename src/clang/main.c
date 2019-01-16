#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{


    const size_t N1 = 20, N2 = 20;
    const double *restrict x1 = linspace(0, 1, N1);
    const double *restrict x2 = linspace(0, 1, N2);

    const double h1 = 1.0/N1;
    const double h2 = 1.0/N2;
    const double eps = 1e-5;

    free(x1);
    free(x2);
    return 0;
}
