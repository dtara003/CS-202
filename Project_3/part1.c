#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) ( BLOCK_LOW((id)+1,p,n)-1 ) 
#define BLOCK_SIZE(id,p,n) (BLOCK_LOW( (id)+1, p, n) - BLOCK_LOW( (id), p, n  ) )
#define BLOCK_OWNER(index,p,n) ( ( ((p)*(index)+1)-1 ) / (n) )
#define MIN(a,b) ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
	int index, id, p, count = 0;
    	double elapsed_time = 0.0;
    	unsigned long long int i, n, low_value, high_value, size, proc0_size, prime, first, global_count = 0;
	char *marked;
    
   	MPI_Init(&argc, &argv);
   	MPI_Barrier(MPI_COMM_WORLD);
   	elapsed_time = -MPI_Wtime();
   	MPI_Comm_rank (MPI_COMM_WORLD, &id);
   	MPI_Comm_size (MPI_COMM_WORLD, &p);
    	
	if (argc != 2) {
          	if (!id) printf ("Command line: %s <m>\n", argv[0]);
          	MPI_Finalize(); exit (1);
    	}
   	
	n = atoll(argv[1]);
   	
	// define start and end values, size
	low_value = 3 + BLOCK_LOW(id,p,n-2) + BLOCK_LOW(id,p,n-2) % 2;
   	high_value = 3 + BLOCK_HIGH(id,p,n-2) - BLOCK_HIGH(id,p,n-2) % 2;
   	size = (high_value - low_value) / 2 + 1;
   	proc0_size = (n-2)/(2*p);
   	marked=(char*)malloc(size);

	// check validity
	if ((3 + proc0_size) < (int) sqrt((double) n)) {
      		if (!id) printf ("Too many processes\n");
      		MPI_Finalize();
      		exit (1);
   	}

   	if (marked == NULL) {
      		printf ("Cannot allocate enough memory\n");
      		MPI_Finalize();
      		exit (1);
   	}

	// set all to prime as default
   	for (i = 0; i < size; i++) marked[i] = 0;
   
	// start at 2 because we already know about all evens
	index = 0;
	prime = 3;
	do {
      		if (prime * prime > low_value) first = (prime * prime - low_value) / 2;
		else {
         		if (!(low_value % prime)) first = 0;
         		else {
				if ((low_value % prime) % 2 == 0) first = prime - (low_value % prime) / 2;
				else first = (prime - (low_value % prime)) / 2;
			}
      		}
      		
		// mark all composites in steps
		for (i = first; i < size; i += prime) marked[i] = 1;
      	
		if (!id) {
         		while (marked[++index]);
         		prime = 2 * index + 3;
     		}
      
		if (p > 1) MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
   	} while (prime * prime <= n);
   	
	// count primes
	for (i = 0; i < size; i++) {
      		if (!marked[i]) count++;
   	}

	if (p > 1) MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	elapsed_time += MPI_Wtime();
   
	if (!id) {
		global_count++;
      		printf ("The total number of primes: %llu, total time: %10.6f", global_count, elapsed_time);
   	}
   	
	MPI_Finalize();
   	
	return 0;
}
