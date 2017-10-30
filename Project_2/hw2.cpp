#include <iostream>
#include "lapacke.h"
#include "blas.h"
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

double randVal() {
    double div = RAND_MAX / 9.0;
    return (rand() / div);
}

void mydgetrf(double* A, double* B, int* pvt, int n) {
    // check contents
    /*for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << A[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    cout << endl;*/
    
    for (int i = 0; i < n - 1; i++) {
        // pivot
        int maxind = i;
        double max = abs(A[i * n + i]);
        
        for (int t = i + 1; t < n; t++) {
            if (abs(A[t * n + i]) > max) {
                maxind = t;
                max = abs(A[t * n + i]);
            }
        }
        
        if (max == 0.0) {
            cout << "LU factorization failed: coefficient matrix is singular" << endl;
            return;
        } else {
            if (maxind != i) {
                // save pivoting information
                int temps = pvt[i];
                pvt[i] = pvt[maxind];
                pvt[maxind] = temps;
                
                double* tempv = (double*)malloc(n * sizeof(double));
                tempv[0] = 0.0;
                
                // swap rows
                for (int j = 0; j < n; j++) {
                    tempv[j] = A[i * n + j];
                    A[i * n + j] = A[maxind * n + j];
                    A[maxind * n + j] = tempv[j];
                }
                
                free(tempv);
            }
        }
        
        // factorization
        for (int j = i + 1; j < n; j++) {
            A[j * n + i] = A[j * n + i] / A[i * n + i];
            
            for (int k = i + 1; k < n; k++) {
                A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
            }
        }
    }
    
    // check contents
    /*cout << "MY LU DECOMP:" << endl << endl << "A = " << endl;
    for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << A[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    cout << endl << "pivot = " << endl;
    for (int i = 0; i < n; i++) {
        cout << pvt[i] << endl;
    }
    cout << endl;*/
    
    return;
}

void mydtrsm1(double* A, double* B, int* pvt, double* y, int n) {
    y[0] = B[pvt[0]];
    
    for (int i = 1; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += y[j] * A[i * n + j];
        }
        y[i] = B[pvt[i]] - sum;
    }
    
    return;
}

void mydtrsm2(double* A, double* x, double* y, int n) {
    x[n - 1] = y[n - 1] / A[(n - 1) * n + (n - 1)];
    
    for (int i = n - 2; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += x[j] * A[i * n + j];
        }
        x[i] = (y[i] - sum) / A[i * n + i];
    }
    
    return;
}

int main()
{
    srand(time(NULL));
	
	int nSize[5] = {1000, 2000, 3000, 4000, 5000};
	
	int x;
    for (x = 0; x < 5; x++) {
        int n = nSize[x];
        char    TRANS = 'N';
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
        int     N = n;
        int     NRHS = 1;
        int     IPIV[n];
        
        cout << "n = " << n << endl << endl;
        
        // allocate arrays
        double* A = (double*)malloc(n * n * sizeof(double));
        double* B = (double*)malloc(n * sizeof(double));
        
        double* A2 = (double*)malloc(n * n * sizeof(double));
        double* B2 = (double*)malloc(n * sizeof(double));
        int* pvt = (int*)malloc(n * sizeof(int));
        double* x = (double*)malloc(n * sizeof(double));
        double* y = (double*)malloc(n * sizeof(double));
        
        // generate random matrices
        for (int i = 0; i < (n * n); i++) {
            if (i < n) {
                A2[i] = randVal();
                B[i] = randVal();
                B2[i] = B[i];
                pvt[i] = i;
                IPIV[i] = i;
            } else {
                A2[i] = randVal();
            }
        }
        // transpose A
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[j * n + i] = A2[i * n + j];
            }
        }
        x[0] = 0.0; y[0] = 0.0;
        
        // output for debugging
        /*cout << "A TRANSPOSED = " << endl;
    	for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << A[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    	cout << endl;
    	cout << "A = " << endl;
    	for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << A2[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    	cout << endl;
    	cout << "B = " << endl;
    	for (int i = 0; i < n; i++) {
    		cout << B[i] << endl;
    	}
    	cout << endl;*/
	
	    clock_t start, end;
	    
	    start = clock();
        // LU factorization
        LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);
        end = clock();
        
        cout << "LAPACK execution time: " << (end - start) / CLOCKS_PER_SEC << endl;
        cout << "GFLOPS: " << (2.0/3) * n * n * n * CLOCKS_PER_SEC / ((end - start) * 1000000000) << endl;
        // check contents
        /*cout << "LAPACK LU DECOMP:" << endl << endl << "A = " << endl;
        for (int i = 0; i < n; i++) {
        		for (int j = 0; j < n; j++) {
        			cout << A[i * n + j] << " ";
        		}
        		cout << endl;
        	}
        cout << endl << "pivot = " << endl;
        for (int i = 0; i < n; i++) {
            cout << IPIV[i] << endl;
        }
        cout << endl;*/
    
        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;
        // This function solve the Ax=B directly
        //dgetrs_(&TRANS,&N,&NRHS,A,&LDA,IPIV,B,&LDB,&INFO);
    
        // change the order of B according to IPIV[] from LU factorization
    
        for(int i = 0; i < N; i++)
        {
            double tmp = B[IPIV[i]-1];
    	B[IPIV[i]-1] = B[i];
    	B[i] = tmp;
        }
        
        // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
        /*cout << "LAPACK FORWARD RESULT" << endl;
	for( int j = 0; j < n; j++) {
		cout << B[j] << " ";
	}
	cout << endl;*/

	    UPLO = 'U';
        DIAG = 'N';
        // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
        
        /*cout << "print the result : { ";
        int i;
        for (i=0;i<N;i++)
        {
    	    cout << B[i] << " ";
        }
        cout << "}" << endl;*/
        
        // LU decomposition
        start = clock();
	mydgetrf(A2, B2, pvt, n);
        end = clock();
	
	cout << "MYDGETRF execution time: " << (end - start) / CLOCKS_PER_SEC << endl;
	cout << "GFLOPS: " << (2.0/3) * n * n * n * CLOCKS_PER_SEC / ((end - start) * 1000000000) << endl; 
        // forward substitution
        mydtrsm1(A2, B2, pvt, y, n);
        /*cout << "MY FORWARD RESULT" << endl;
    	for (int j = 0; j < n; j++) {
    		cout << y[j] << " ";
    	}
    	cout << endl;*/
    	// backward substitution
        mydtrsm2(A2, x, y, n);
        
        /*cout << endl;
        for (int i = 0; i < n; i++) {
            cout << x[i] << endl;
        }*/
        
	double maxDiff = 0.0;
	for (int i = 0; i < n; i++) {
		double tempDiff = abs(B[i] - x[i]);
		if (tempDiff > maxDiff) {
			maxDiff = tempDiff;
		}
	}
	cout << "Error check - max difference: " << maxDiff << endl;
	if (maxDiff < 0.0001) {
		cout << "Error < 1e-3. Negligable." << endl;
	} else {
		cout << "ERROR GREATAER THAN 1e-3." << endl;
	}
	cout << endl;
	
	free(A); free(B); free(A2); free(B2); free(pvt); free(x); free(y);
    }
    
    return 0;
}
