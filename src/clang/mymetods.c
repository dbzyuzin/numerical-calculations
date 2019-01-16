
double* linspace(int xa, int xb, int N) {
	double *x = calloc(N, sizeof(double));
    double h = (xb - xa)/(double)(N-1);
    for (int i=0; i<N; i++){
        x[i] = xa + i*h;
    }
	return x;
}

void print_res(int N1,int N2,double h1,doube h2,double eps,int iter_count,double rka) {
    printf("Параметры :\n%7s %7s %7s %7s %8s\n","N1", "N2", "h1", "h2", "eps");
    // printf("\n%7d %7d %7.3f %7.3f %8e\n","N1", "N2", "h1", "h2", "eps");
    // s += "{: ^7d}{: ^7d}{: ^7.3f}{: ^7.3f}{: ^8}\n\n".format(N1, N2, h1, h2, eps)
    // s += "Результаты:\n {: ^12}{: ^24}\n".format("Iter count", "Max Fail")
    // s += " {: ^12}{: ^24.15f}\n{}\n".format(iter_count, np.max(yerr), "="*36)
}
