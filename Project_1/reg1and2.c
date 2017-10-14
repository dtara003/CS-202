#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double randVal() {
    double div = RAND_MAX / 9.0;
    return (rand() / div);
}

int main() {
    srand(time(NULL));
    // printf("%f", randVal());
    
    int nSize[6] = {64, 128, 256, 512, 1024, 2048};
    
    int x;
    for (x = 0; x < 6; x++) {
        int n = nSize[x];
        printf("n = %d", n);
        printf("\n\n");
        
        int i, j, k;
        
        // allocate arrays
        double** a = (double**)malloc(n * sizeof(double*));
        for (i = 0; i < n; i++) {
            a[i] = (double*)malloc(n * sizeof(double));
        }
        double** b = (double**)malloc(n * sizeof(double*));
        for (i = 0; i < n; i++) {
            b[i] = (double*)malloc(n * sizeof(double));
        }
        double** c1 = (double**)malloc(n * sizeof(double*));
        for (i = 0; i < n; i++) {
            c1[i] = (double*)malloc(n * sizeof(double));
        }
        double** c2 = (double**)malloc(n * sizeof(double*));
        for (i = 0; i < n; i++) {
            c2[i] = (double*)malloc(n * sizeof(double));
        }
        double** c3 = (double**)malloc(n * sizeof(double*));
        for (i = 0; i < n; i++) {
            c3[i] = (double*)malloc(n * sizeof(double));
        }
        
        // generate random matrices
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                a[i][j] = randVal();
                b[i][j] = randVal();
                double temp = randVal();
                c1[i][j] = temp;
                c2[i][j] = temp;
                c3[i][j] = temp;
            }
        }
        
        clock_t start, end;
        
        start = clock();
        
        // dgemm0: simple ijk triple loop algorithm
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    c1[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        
        end = clock();
        
        printf("dgemm0 in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
        printf("\n");
        
        start = clock();
        
        // dgemm1: simple ijk triple loop algorithm with register reuse
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                register double r = c2[i][j];
                
                for (k = 0; k < n; k++) {
                    r += a[i][k] * b[k][j];
                }
                
                c2[i][j] = r;
            }
        }
        
        end = clock();
        
        printf("dgemm1 in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
        printf("\n");
        
        start = clock();
        
        // dgemm2: using 12 registers
        for (i = 0; i < n; i += 2) {
            for (j = 0; j < n; j += 2) {
                register double cr1 = c3[i][j];
                register double cr2 = c3[i + 1][j];
                register double cr3 = c3[i][j + 1];
                register double cr4 = c3[i + 1][j + 1];
                
                for (k = 0; k < n; k += 2) {
                    //registers for a[...] and b[...]
                    register double a1 = a[i][k];
                    register double a2 = a[i + 1][k];
                    register double a3 = a[i][k + 1];
                    register double a4 = a[i + 1][k + 1];
                    
                    register double b1 = b[k][j];
                    register double b2 = b[k + 1][j];
                    register double b3 = b[k][j + 1];
                    register double b4 = b[k + 1][j + 1];
                    
                    cr1 += a1 * b1 + a3 * b2;
                    cr2 += a2 * b1 + a4 * b2;
                    cr3 += a1 * b3 + a3 * b4;
                    cr4 += a2 * b3 + a4 * b4;
                }
                
                c3[i][j] = cr1;
                c3[i + 1][j] = cr2;
                c3[i][j + 1] = cr3;
                c3[i + 1][j + 1] = cr4;
            }
        }
        
        end = clock();
        
        printf("dgemm2 in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
        printf("\n\n");
        
        // max difference between c1 and c2 and between c1 and c3
        double diff1 = 0;
        double diff2 = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double temp = abs(c1[i][j] - c2[i][j]);
                if (temp > diff1) {
                    diff1 = temp;
                }
                
                temp = abs(c1[i][j] - c3[i][j]);
                if (temp > diff2) {
                    diff2 = temp;
                }
            }
        }
        
        printf("MAX DIFFERENCE:\ncomparing dgemm0 and dgemm1: %f", diff1);
        printf("\ncomparing dgemm0 and dgemm2: %f \n\n", diff2);
        
        // free memory
        for (i = 0; i < n; i++) {
            free(a[i]);
            free(b[i]);
            free(c1[i]);
            free(c2[i]);
            free(c3[i]);
        }
        free(a); free(b); free(c1); free(c2); free(c3);
    }
    
    return 0;
}
