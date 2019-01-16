
void linspace(double *x, int xa, int xb, int N) {
    double h = (xb - xa)/(double)(N-1);
    for (int i=0; i<N; i++){
        x[i] = xa + i*h;
    }
}
