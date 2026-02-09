/*#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(int argc, char *argv[]) {

    int threads = 1;

    if (argc > 1) {
        threads = atoi(argv[1]);
    }

    omp_set_num_threads(threads);

    int N = 1000;
    double epsilon = 1.0, sigma = 1.0;
    double rc2 = 2.5 * 2.5;

    double *x = malloc(N * sizeof(double));
    double *y = malloc(N * sizeof(double));
    double *z = malloc(N * sizeof(double));

    for (int i = 0; i < N; i++) {
        x[i] = drand48();
        y[i] = drand48();
        z[i] = drand48();
    }

    double potential = 0.0;
    double start = omp_get_wtime();

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];

            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-4 || r2 > rc2) continue;

            double inv2 = (sigma*sigma) / r2;
            double inv6 = inv2 * inv2 * inv2;
            double inv12 = inv6 * inv6;

            potential += 4 * epsilon * (inv12 - inv6);
        }
    }

    double end = omp_get_wtime();

    printf("Threads Used: %d\n", threads);
    printf("Execution Time: %f seconds\n", end - start);
    printf("Total Potential Energy: %f\n", potential);

    free(x); free(y); free(z);
    return 0;
}

*/


//Q2

/*
#include <stdio.h>
#include <string.h>
#include <omp.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

int main() {

    char A[] = "ACACACTA";
    char B[] = "AGCACACA";

    int m = strlen(A);
    int n = strlen(B);

    int H[m+1][n+1];

    int match = 2;
    int mismatch = -1;
    int gap = -1;


    for (int i = 0; i <= m; i++)
        for (int j = 0; j <= n; j++)
            H[i][j] = 0;

    double start = omp_get_wtime();


    for (int d = 1; d <= m + n - 1; d++) {

        #pragma omp parallel for
        for (int i = 1; i <= m; i++) {

            int j = d - i;

            if (j >= 1 && j <= n) {

                int score;
                if (A[i-1] == B[j-1])
                    score = match;
                else
                    score = mismatch;

                int diag = H[i-1][j-1] + score;
                int up   = H[i-1][j] + gap;
                int left = H[i][j-1] + gap;

                H[i][j] = MAX(0, MAX(diag, MAX(up, left)));
            }
        }
    }

    double end = omp_get_wtime();

    printf("Smith-Waterman completed\n");
    printf("Execution Time: %f seconds\n", end - start);

    return 0;
}


*/


#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char *argv[]) {

    int N = 1500;
    int steps = 1500;
    int threads = 1;

    if (argc > 1)
        threads = atoi(argv[1]);

    omp_set_num_threads(threads);

    double **T = malloc(N * sizeof(double*));
    double **Tnew = malloc(N * sizeof(double*));

    for (int i = 0; i < N; i++) {
        T[i] = calloc(N, sizeof(double));
        Tnew[i] = calloc(N, sizeof(double));
    }

    for (int j = 0; j < N; j++)
        T[0][j] = 100.0;

    double start = omp_get_wtime();

    for (int t = 0; t < steps; t++) {
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) {
            for (int j = 1; j < N-1; j++) {
                Tnew[i][j] = 0.25 * (
                    T[i+1][j] + T[i-1][j] +
                    T[i][j+1] + T[i][j-1]
                );
            }
        }

        double **temp = T;
        T = Tnew;
        Tnew = temp;
    }

    double end = omp_get_wtime();

    printf("Threads Used: %d\n", threads);
    printf("Execution Time: %f seconds\n", end - start);

    for (int i = 0; i < N; i++) {
        free(T[i]);
        free(Tnew[i]);
    }
    free(T);
    free(Tnew);

    return 0;
}

