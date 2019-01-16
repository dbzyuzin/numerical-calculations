#include <stdlib.h>
#include "mymetods.h"

void edge_computing(double *ys, double (*u)(const double, const double), 
            const double* restrict x1, const size_t N1, 
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