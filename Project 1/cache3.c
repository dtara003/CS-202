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
    
    int n = 64;
    printf("n = %d", n);
    printf("\n");
    int B = 32;
    printf("B = %d", B);
    printf("\n\n");
    
    int i, j, k, i1, j1, k1;
    
    double error[12];
    error[0] = 0.0;
    
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
    double** cSTART = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        cSTART[i] = (double*)malloc(n * sizeof(double));
    }
    
    // generate random matrices
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            a[i][j] = randVal();
            b[i][j] = randVal();
            
            double temp = randVal();
            cSTART[i][j] = temp;
            c1[i][j] = temp;
            c2[i][j] = temp;
            c3[i][j] = temp;
        }
    }
    
    clock_t start, end;
    
    printf("SIMPLE TRIPLE LOOP W/ REG REUSE\n\n");
    
    start = clock();
    
    // simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            register double r = c3[i][j];
            
            for (k = 0; k < n; k++) {
                r += a[i][k] * b[k][j];
            }
            
            c3[i][j] = r; // used to check error
        }
    }
    
    end = clock();
    
    printf("IJK in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // clear c1 for use again
    /*for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c1[i][j] = 0;
        }
    }*/
    
    start = clock();
    
    // simple JIK
    for (j = 0; j < n; j++) {
        for (i = 0; i < n; i++) {
            register double r = c1[i][j];
            
            for (k = 0; k < n; k++) {
                r += a[i][k] * b[k][j];
            }
            
            c1[i][j] = r;
        }
    }
    
    end = clock();
    
    printf("JIK in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c1[i][j]);
            if (temp > error[1]) {
                error[1] = temp;
            }
        }
    }
    
    // clear c1 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c1[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // simple KIJ
    for (k = 0; k < n; k++) {
        for (i = 0; i < n; i++) {
            register double r = a[i][k];
            
            for (j = 0; j < n; j++) {
                c1[i][j] += r * b[k][j];
            }
        }
    }
    
    end = clock();
    
    printf("KIJ in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c1[i][j]);
            if (temp > error[2]) {
                error[2] = temp;
            }
        }
    }
    
    // clear c1 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c1[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // simple IKJ
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            register double r = a[i][k];
            
            for (j = 0; j < n; j++) {
                c1[i][j] += r * b[k][j];
            }
        }
    }
    
    end = clock();
    
    printf("IKJ in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c1[i][j]);
            if (temp > error[3]) {
                error[3] = temp;
            }
        }
    }
    
    // clear c1 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c1[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // simple JKI
    for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++) {
            register double r = b[k][j];
            
            for (i = 0; i < n; i++) {
                c1[i][j] += a[i][k] * r;
            }
        }
    }
    
    end = clock();
    
    printf("JKI in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c1[i][j]);
            if (temp > error[4]) {
                error[4] = temp;
            }
        }
    }
    
    // clear c1 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c1[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // simple KJI
    for (k = 0; k < n; k++) {
        for (j = 0; j < n; j++) {
            register double r = b[k][j];
            
            for (i = 0; i < n; i++) {
                c1[i][j] += a[i][k] * r;
            }
        }
    }
    
    end = clock();
    
    printf("KJI in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c1[i][j]);
            if (temp > error[5]) {
                error[5] = temp;
            }
        }
    }
    
    printf("BLOCKED VERSIONS\n\n");
    
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
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c2[i][j]);
            if (temp > error[6]) {
                error[6] = temp;
            }
        }
    }
    
    // clear c2 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c2[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // blocked JIK
    for (j = 0; j < n; j += B) {
        for (i = 0; i < n; i += B) {
            for (k = 0; k < n; k += B) {
                // B x B mini matrix multiplications
                for (j1 = j; j1 < j + B; j1++) {
                    for (i1 = i; i1 < i + B; i1++) {
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
    
    printf("blocked JIK in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c2[i][j]);
            if (temp > error[7]) {
                error[7] = temp;
            }
        }
    }
    
    // clear c2 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c2[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // blocked KIJ
    for (k = 0; k < n; k += B) {
        for (i = 0; i < n; i += B) {
            for (j = 0; j < n; j += B) {
                // B x B mini matrix multiplications
                for (k1 = k; k1 < k + B; k1++) {
                    for (i1 = i; i1 < i + B; i1++) {
                        register double r = a[i1][k1];
                        
                        for (j1 = j; j1 < j + B; j1++) {
                            c2[i1][j1] += r * b[k1][j1];
                        }
                    }
                }
            }
        }
    }
    
    end = clock();
    
    printf("blocked KIJ in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c2[i][j]);
            if (temp > error[8]) {
                error[8] = temp;
            }
        }
    }
    
    // clear c2 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c2[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // blocked KIJ
    for (i = 0; i < n; i += B) {
        for (k = 0; k < n; k += B) {
            for (j = 0; j < n; j += B) {
                // B x B mini matrix multiplications
                for (i1 = i; i1 < i + B; i1++) {
                    for (k1 = k; k1 < k + B; k1++) {
                        register double r = a[i1][k1];
                        
                        for (j1 = j; j1 < j + B; j1++) {
                            c2[i1][j1] += r * b[k1][j1];
                        }
                    }
                }
            }
        }
    }
    
    end = clock();
    
    printf("blocked IKJ in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c2[i][j]);
            if (temp > error[9]) {
                error[9] = temp;
            }
        }
    }
    
    // clear c2 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c2[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // blocked JKI
    for (j = 0; j < n; j += B) {
        for (k = 0; k < n; k += B) {
            for (i = 0; i < n; i += B) {
                // B x B mini matrix multiplications
                for (j1 = j; j1 < j + B; j1++) {
                    for (k1 = k; k1 < k + B; k1++) {
                        register double r = b[k1][j1];
                        
                        for (i1 = i; i1 < i + B; i1++) {
                            c2[i1][j1] += a[i1][k1] * r;
                        }
                    }
                }
            }
        }
    }
    
    end = clock();
    
    printf("blocked JKI in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c2[i][j]);
            if (temp > error[10]) {
                error[10] = temp;
            }
        }
    }
    
    // clear c2 for use again
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c2[i][j] = cSTART[i][j];
        }
    }
    
    start = clock();
    
    // blocked KJI
    for (k = 0; k < n; k += B) {
        for (j = 0; j < n; j += B) {
            for (i = 0; i < n; i += B) {
                // B x B mini matrix multiplications
                for (k1 = k; k1 < k + B; k1++) {
                    for (j1 = j; j1 < j + B; j1++) {
                        register double r = b[k1][j1];
                        
                        for (i1 = i; i1 < i + B; i1++) {
                            c2[i1][j1] += a[i1][k1] * r;
                        }
                    }
                }
            }
        }
    }
    
    end = clock();
    
    printf("blocked KJI in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
    printf("\n\n");
    
    // check error against simple IJK
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c3[i][j] - c2[i][j]);
            if (temp > error[11]) {
                error[11] = temp;
            }
        }
    }
    
    // max difference between c1 and c2 and between c1 and c3
    double diff1 = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double temp = abs(c1[i][j] - c2[i][j]);
            if (temp > diff1) {
                diff1 = temp;
            }
        }
    }
    
    printf("MAX DIFFERENCE (compared to simple IJK implementation:\n\n");
    printf("simple IJK: %f\n", error[0]); printf("simple JIK: %f\n", error[1]);
    printf("simple KIJ: %f\n", error[2]); printf("simple IKJ: %f\n", error[3]);
    printf("simple JKI: %f\n", error[4]); printf("simple KJI: %f\n", error[5]);
    printf("simple IJK: %f\n", error[6]); printf("simple JIK: %f\n", error[7]);
    printf("simple KIJ: %f\n", error[8]); printf("simple IKJ: %f\n", error[9]);
    printf("simple JKI: %f\n", error[10]); printf("simple KJI: %f\n", error[11]);
    
    // output C1, C2, C3 for debugging
    /*for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", c1[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", c2[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", c3[i][j]);
        }
        printf("\n");
    }
    
    printf("\n");*/
    
    // free memory
    for (i = 0; i < n; i++) {
        free(a[i]);
        free(b[i]);
        free(c1[i]);
        free(c2[i]);
        free(c3[i]);
        free(cSTART[i]);
    }
    
    free(a); free(b); free(c1); free(c2); free(c3); free(cSTART);
    
    return 0;
}
