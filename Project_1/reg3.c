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
    int mSize[6] = {66, 129, 258, 513, 1026, 2049};
    
    int x;
    for (x = 0; x < 6; x++) {
        int n = nSize[x];
        printf("n = %d", n);
        printf("\n\n");
        int m = mSize[x];
        
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
        double** c3 = (double**)malloc(n * sizeof(double*));
        for (i = 0; i < n; i++) {
            c3[i] = (double*)malloc(n * sizeof(double));
        }
        
        double** aa = (double**)malloc(m * sizeof(double*));
        for (i = 0; i < m; i++) {
            aa[i] = (double*)malloc(m * sizeof(double));
        }
        double** bb = (double**)malloc(m * sizeof(double*));
        for (i = 0; i < m; i++) {
            bb[i] = (double*)malloc(m * sizeof(double));
        }
        double** c4 = (double**)malloc(m * sizeof(double*));
        for (i = 0; i < m; i++) {
            c4[i] = (double*)malloc(m * sizeof(double));
        }
        
        c3[0][0] = 0.0;
        aa[0][0] = 0.0;
        bb[0][0] = 0.0;
        c4[0][0] = 0.0;
        
        // generate random matrices
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                a[i][j] = randVal();
                b[i][j] = randVal();
                
                aa[i][j] = a[i][j];
                bb[i][j] = b[i][j];
            }
        }
        
        // print contents of A and AA
        /*for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%f ", a[i][j]);
            }
            printf("\n");
        }
        printf("\n");
        for (i = 0; i < m; i++) {
            for (j = 0; j < m; j++) {
                printf("%f ", aa[i][j]);
            }
            printf("\n");
        }
        printf("\n");*/
        
        clock_t start, end;
        
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
        
        start = clock();
        
        // dgemm3: max of 16 registers
        /*for (i = 0; i < n; i += 2) {
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
        }*/
        
        end = clock();
        
        printf("dgemm3 in seconds: %f", ((double)(end - start)) / CLOCKS_PER_SEC);
        printf("\n\n");
        
        // output C
        /*printf("\n");
        
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%f ", c3[i][j]);
            }
            printf("\n");
        }
        
        printf("\n");
        
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%f ", c4[i][j]);
            }
            printf("\n");
        }
        
        printf("\n");*/
        
        // max difference between c3 and c4
        /*double diff1 = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double temp = abs(c3[i][j] - c4[i][j]);
                if (temp > diff1) {
                    diff1 = temp;
                }
            }
        }
        
        printf("max error:\ncomparing dgemm0 and dgemm1: %f", diff1);
        printf("\ncomparing dgemm0 and dgemm2: %f", diff2);
        
        if (diff1 < 0.001) {
            printf("\nCORRECT\n\n");
        } else {
            printf("\nNOT CORRECT\n\n");
        }*/
        
        // free memory
        for (i = 0; i < n; i++) {
            free(a[i]);
            free(b[i]);
            free(c3[i]);
        }
        free(a); free(b); free(c3);
        
        for (i = 0; i < m; i++) {
            free(aa[i]);
            free(bb[i]);
            free(c4[i]);
        }
        free(aa); free(bb); free(c4);
    }
    
    return 0;
}
