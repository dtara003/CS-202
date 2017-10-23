#include <iostream>
#include "lapacke.h"
#include "blas.h"
using namespace std;

double randVal() {
    double div = RAND_MAX / 9.0;
    return (rand() / div);
}

int main()
{
    srand(time(NULL));
    
    /*// use new to allocate memory if you need large space
    // Here, we want to solve AX = b
    //    x1 + 2x2 + 3x3 = 1
    //    2x1 + x2 + x3  = 1
    //    x1 + x2 + x3   = 1
    // in C, you should initialize A as:
    //  A = { 1 2 3
    //        2 1 1
    //        1 1 1 }
    // IF you use this A to call LAPACK function, it gets a wrong result
    
    // BUT, LAPACK need the A to store in COLUMN-order
    // SO, we initial A as (for the same system):
    //  A' = { 1 2 1
    //         2 1 1
    //         3 1 1 }
    // correct solution = {0 2 -1}'*/
    
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
        double* A = (double*)malloc(n * n * sizeof(double));
        double* B = (double*)malloc(n * sizeof(double));
        
        // generate random matrices
        for (int i = 0; i < (n * n); i++) {
            if (i < n) {
                A[i] = randVal();
                B[i] = randVal();
            } else {
                A[i] = randVal();
            }
        }
        
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
        
        cout << "print the result : {";
        int i;
        for (i=0;i<N;i++)
        {
    	cout << B[i] << " ";
        }
        cout << "}" << endl;
    }
    
    return 0;
}
