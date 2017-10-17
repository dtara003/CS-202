#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double randVal() {
    double div = RAND_MAX / 9.0;
    return (rand() / div);
}

int main() {
    srand(time(NULL));
    
    int n = 2048;
    printf("n = %d", n);
    printf("\n");
    int B = 16;
    printf("B = %d", B);
    printf("\n\n");
    
    int i, j, k, i1, j1, k1;
    
    // allocate arrays
    double** a = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        a[i] = (double*)malloc(n * sizeof(double));
    }
    double** b = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        b[i] = (double*)malloc(n * sizeof(double));
    }
    double** c2 = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        c2[i] = (double*)malloc(n * sizeof(double));
    }
    double** c3 = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        c3[i] = (double*)malloc(n * sizeof(double));
    }
    
    clock_t start, end;
    
    start = clock();
    
    // blocked IJK 
    for (i = 0; i < n; i += B) {
        for (j = 0; j < n; j += B) {
            for (k = 0; k < n; k += B) {
                // B x B mini matrix multiplications
                for (i1 = i; i1 < i + B; i1++) {
                    for (j1 = j; j1 < j + B; j1++) {
                        register double r = c2[i1][j1];
                        
                        for (k1 = k; k1 < k + B; k1++) {
                            r += a[i1][k1] * b[k1][j1];
                        }
                        
                        c2[i1][j1] = r;
                    }    
                }
            }
        }
    }
    
    end = clock();
    
    printf("blocked IJK in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    printf("GFLOPS: %.16f", 2*((double)pow(n,3))/((((double)(end - start)))*pow(10,9)));
    printf("\n");
    
    
    start = clock();
        
    // reg and cache reuse
    for (i = 0; i < n; i += B) {
        for (j = 0; j < n; j += B) {
            for (k = 0; k < n; k += B) {
                for (i1 = i; i1 < i + B; i1 += 2) {
                    for (j1 = j; j1 < j + B; j1 += 2) {
                        register double cr1 = c3[i][j];
                        register double cr2 = c3[i + 1][j];
                        register double cr3 = c3[i][j + 1];
                        register double cr4 = c3[i + 1][j + 1];
                        
                        for (k1 = k; k1 < k + B; k1 += 2) {
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
            }
        }
    }
        
    end = clock();
    
    printf("combined cache and register reuse for IJK in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    printf("GFLOPS: %.16f", 2*((double)pow(n,3))/((((double)(end - start)))*pow(10,9)));
    printf("\n");
    
    // max difference between c1 and c2 and between c1 and c3
    double diff1 = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c2[i][j] - c3[i][j]);
            if (temp > diff1) {
                diff1 = temp;
            }
        }
    }
    printf("MAX DIFFERENCE: %f\n", diff1);
    
    // free memory
    for (i = 0; i < n; i++) {
        free(a[i]);
        free(b[i]);
        free(c2[i]);
        free(c3[i]);
    }
    
    free(a); free(b); free(c2); free(c3);
    
    return 0;
}
