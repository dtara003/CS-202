#include <iostream>
#include "lapacke.h"
#include "blas.h"
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

// generates a random value
double randVal() {
    double div = RAND_MAX / 9.0;
    return (rand() / div);
}

// part 1.2 implementation of LU decomposition and merge
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
void dgemm3(double* A, double* B, double* C, int i, int k, int n, int block) {
    int j = i;
    int k2 = i;
    
    // modified IJK BLOCK algorithm
    for (i; i < n; i += block) {
        for (k; j < j2; j += block) {
            for (k; k < n; k += block) {
                for (int m = i; m < i + block && m < n; m += 4) {
                    for (int n = k; n < k + block && n < n; n += 2) {
                        // pull C to register and then store back into C later
                        register double c0 = C[m * n + n];               register double c4 = C[(m + 2) * n + n];
                        register double c1 = C[m * n + (n + 1)];         register double c5 = C[(m + 2) * n + (n + 1)];
                        register double c2 = C[(m + 1) * n + n];         register double c6 = C[(m + 3) * n + n];
                        register double c3 = C[(m + 1) * n + (n + 1)];   register double c7 = C[(m + 3) * n + (n + 1)];
                        
                        for (int p = j; p < j + block && p < k2; p += 2) {
                            // pull A and B portions to registers
                            // performing Ijj on 4 x 2 matrix
                            register double a0 = A[m * n + p];            register double a2 = A[(m + 2) * n + p];
                            register double a1 = A[(m + 1) * n + p];      register double a3 = A[(m + 3) * n + p];
                            register double b0 = B[p * n + n];            register double b1 = B[p * n + (n + 1)];
                            
                            // subtract from previous C values
                            c0 = c0 - a0 * b0;    c1 = c1 - a0 * b1;
                            // new A value because old one done being used
                            a0 = A[m * n + (p + 1)];
                            c2 = c2 - a1 * b0;    c3 = c3 - a1 * b1;
                            // new A value because old one done being used
                            a1 = A[(m + 1) * n + p + 1];
                            c4 = c4 - a2 * b0;    c5 = c5 - a2 * b1;
                            // new A value because old one done being used
                            a2 = A[(m + 2) * n + p + 1];
                            c6 = c4 - a3 * b0;    c7 = c7 - a3 * b1;
                            // new A value because old one done being used
                            a3 = A[(m + 3) * n + p + 1];
                            // new B value because old one done being used
                            b0 = B[(p + 1) * n + n]; b1 = B[(p + 1) * n + (n + 1)];
                            
                            c0 = c0 - a0 * b0;    c1 = c1 - a0 * b1;
                            c2 = c2 - a1 * b0;    c3 = c3 - a1 * b1;
                            c4 = c4 - a2 * b0;    c5 = c5 - a2 * b1;
                            c6 = c4 - a3 * b0;    c7 = c7 - a3 * b1;
                        }
                        // store back into C
                        C[m * n + n] = c0;               C[(m + 2) * n + n] = c4;
                        C[m * n + (n + 1)] = c1;         C[(m + 2) * n + (n + 1)] = c5;
                        C[(m + 1) * n + n] = c2;         C[(m + 3) * n + n] = c6;
                        C[(m + 1) * n + (n + 1)] = c3;   C[(m + 3) * n + (n + 1)] = c7;
                    }
                }
            }
        }
    }
    
    /*for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << C[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    cout << endl;*/
    
    return;
}
void myblockeddgetrf(double* A, int* pvt, int n, int block) {
    // check contents
    /*for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << A[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    cout << endl;*/
    
    double max;
    
    for (int i = 0; i < n; i += block) {
        for (int j = i; j < i + block && j < n; j++) {
            int maxind = j;
            max = abs(A[j * n + j]);
            
            for (int t = j + 1; t < n; t++) {
                if (abs(A[t * n + j]) > max) {
                    maxind = t;
                    max = abs(A[t * n + j]);
                }
            }
            
            // check validity before moving on
            if (max == 0.0) {
                cout << "LU factorization failed: coefficient matrix is singular" << endl;
                return;
            } else {
                if (maxind != j) {
                    // store pivot
                    int temps = pvt[j];
                    pvt[j] = pvt[maxind];
                    pvt[maxind] = temps;
                    
                    // perform swap on rows
                    double* tempv = (double*)malloc(n * sizeof(double));
                    tempv[0] = 0.0;
                    
                    for (int k = 0; k < n; k++) {
                        tempv[k] = A[j * n + k];
                        A[j * n + k] = A[maxind * n + k];
                        A[maxind * n + k] = tempv[k];
                    }
                    
                    free(tempv);
                }
            }
            
            for (int k = j + 1; k < n; k++) {
                A[k * n + j] = A[k * n + j] / A[j * n + j];
                
                for (int m = j + 1; m < j + block && m < n; m++) {
                    A[k * n + m] = A[k * n + m] - A[k * n + j] * A[j * n + m];
                }
            }
            
            // check contents
            /*cout << "CURRENT DECOMP:" << endl << endl << "A = " << endl;
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
        }
        int tmp;
        if (i + block < n) {
            tmp = i + block;
        } else {
            tmp = n;
        }
        
        for (int j = i; j < i + block - 1 && j < n - 1; j++) {
            for (int k = j; k < i + block && k < n; k++) {
                for (int m = tmp; m < n; m++) {
                    A[k * n + m] -= A[k * n + j] * A[j * n + m];
                }
                
                // A[i1 * n + j] -= sum;
            }
        }
        
        dgemm3(A, A, A, tmp, i, n, 60);
        
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

int main()
{
    srand(time(NULL));
	
	int nSize[5] = {1000, 2000, 3000, 4000, 5000};
	
	int x;
    for (x = 0; x < 3; x++) {
        int block = 10;
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
        double* A3 = (double*)malloc(n * n * sizeof(double));
        double* B3 = (double*)malloc(n * sizeof(double));
        int* pvt = (int*)malloc(n * sizeof(int));
        double* x = (double*)malloc(n * sizeof(double));
        double* y = (double*)malloc(n * sizeof(double));
        
        // generate random matrices
        for (int i = 0; i < (n * n); i++) {
            if (i < n) {
                A2[i] = randVal();
                A3[i] = A2[i];
                B[i] = randVal();
                B2[i] = B[i];
                B3[i] = B[i];
                pvt[i] = i;
                IPIV[i] = i;
            } else {
                A2[i] = randVal();
                A3[i] = A2[i];
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
        cout << "LAPACK execution time: " << (end - start) / (double)CLOCKS_PER_SEC << endl;
        cout << "GFLOPS: " << (2.0/3) * n * n * n * (double)CLOCKS_PER_SEC / ((end - start) * 1000000000) << endl;
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
        for(int i = 0; i < N; i++) {
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
        
        // PART 1.2 IMPLEMENTATION
        
        // LU decomposition
        start = clock();
	    mydgetrf(A2, B2, pvt, n);
        end = clock();
    	cout << "MYDGETRF execution time: " << (end - start) / (double)CLOCKS_PER_SEC << endl;
    	cout << "GFLOPS: " << (2.0/3) * n * n * n * (double)CLOCKS_PER_SEC / ((end - start) * 1000000000) << endl; 
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
    	if (maxDiff < 0.0001) {
    	    cout << "Error less than 1e-3." << endl;
    	} else {
    		cout << "ERROR GREATER THAN 1e-3." << endl;
    	}
    	
    	// clear pvt and x and y for use again
    	for (int i = 0; i < n; i++) {
    	    pvt[i] = i; x[i] = 0.0; y[i] = 0.0;
    	}
    	
    	start = clock();
    	myblockeddgetrf(A3, pvt, n, block);
    	end = clock();
        cout << "BLOCKED execution time: " << (end - start) / (double)CLOCKS_PER_SEC << endl;
    	cout << "GFLOPS: " << (2.0/3) * n * n * n * (double)CLOCKS_PER_SEC / ((end - start) * 1000000000) << endl; 
        mydtrsm1(A3, B3, pvt, y, n);
        mydtrsm2(A3, x, y, n);
        
        maxDiff = 0.0;
    	for (int i = 0; i < n; i++) {
    		double tempDiff = abs(B[i] - x[i]);
    		if (tempDiff > maxDiff) {
    			maxDiff = tempDiff;
    		}
    	}
    	if (maxDiff < 0.0001) {
    	    cout << "Error less than 1e-3." << endl;
    	} else {
    		cout << "ERROR GREATER THAN 1e-3." << endl;
    	}
    	cout << endl;
    	
    	free(A); free(B); free(A2); free(B2); free(A3); free(B3); free(pvt); free(x); free(y);
    }
    
    return 0;
}