#include <iostream>
#include "lapacke.h"
#include "blas.h"
using namespace std;

double randVal() {
    double div = RAND_MAX / 9.0;
    return (rand() / div);
}

void mydgetrf(double* A, double* B) {
    int n = 2;
    for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << A[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    cout << endl;
    
    return;
}

int main()
{
    srand(time(NULL));
    
    /*double  A[9] =
	{
	    1, 2, 1,
	    2, 1, 1,
	    3, 1, 1
	};
    
    double B[3] =
	{
	    1,
	    1,
	    1
	};*/
	
	int nSize[6] = {2, 1000, 2000, 3000, 4000, 5000};
	
	int x;
    for (x = 0; x < 1; x++) {
        int n = nSize[x];
        char    TRANS = 'N';
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
        int     N = n;
        int     NRHS = 1;
        int     IPIV[n];
        
        cout << "n = " << n << endl;
        
        // allocate arrays
        double* A = (double*)malloc(n * n * sizeof(double));
        double* B = (double*)malloc(n * sizeof(double));
        double* A2 = (double*)malloc(n * n * sizeof(double));
        double* B2 = (double*)malloc(n * sizeof(double));
        
        // generate random matrices
        for (int i = 0; i < (n * n); i++) {
            if (i < n) {
                A[i] = randVal();
                B[i] = randVal();
                A2[i] = A[i];
                B2[i] = B[i];
            } else {
                A[i] = randVal();
                A2[i] = A[i];
            }
        }
        
        // output for debugging
    	for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			cout << A[i * n + j] << " ";
    		}
    		cout << endl;
    	}
    	cout << endl;
    
    	for (int i = 0; i < n; i++) {
    		cout << B[i] << endl;
    	}
	
        // LU factorization
        LAPACK_dgetrf(&N,&N,A,&LDA,IPIV,&INFO);
    
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
        UPLO = 'U';
        DIAG = 'N';
        // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,A, &N, B, &N);
        
        cout << "print the result : { ";
        int i;
        for (i=0;i<N;i++)
        {
    	cout << B[i] << " ";
        }
        cout << "}" << endl;
        
        mydgetrf(A2, B2);
        
        free(A); free(B); free(A2); free(B2);
    }
    
    return 0;
}