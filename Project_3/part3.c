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
	unsigned long long int currHigh, currLow, currSize, currFirst = 0;
	char *currMarked;
	unsigned long long int low3, high3, size3, start3, end3, cacheSize = 0;
	cacheSize = 2000000; start3 = 0;
    
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
	currLow = 3;
	currHigh = (unsigned long long int)sqrt(n);
	currSize = (currHigh - currLow) / 2 + 1;
	currMarked = (char*)malloc(currSize);

	if (currMarked == NULL) {
		printf("Cannot allocate memory\n");
		MPI_Finalize(); exit(1);
	}

	for (i = 0; i < currSize; i++) currMarked[i] = 0;

   	low_value = 3 + BLOCK_LOW(id,p,n-2) + (BLOCK_LOW(id,p,n-2) % 2);
   	high_value = 3 + BLOCK_HIGH(id,p,n-2) - (BLOCK_HIGH(id,p,n-2) % 2);
   	size = (high_value - low_value) / 2 + 1;
   	proc0_size = (n-2)/(2*p);
  
   	if ((3 + proc0_size) < (int) sqrt((double) n)) {
      		if (!id) printf ("Too many processes\n");
      		MPI_Finalize();
      		exit (1);
   	}
   	
	marked = (char *) malloc (size);
   	if (marked == NULL) {
      		printf ("Cannot allocate enough memory\n");
      		MPI_Finalize();
      		exit (1);
   	}
   
	for (i = 0; i < size; i++) marked[i] = 0;
   
	index = 0;
   	prime = 3;

	do {
		currFirst = prime * prime - currLow;
		currFirst /= 2;
		
		for (i = currFirst; i < currSize; i += prime) currMarked[i] = 1;

		while (currMarked[++index]);
		prime = 2 * index + 3;
	} while (prime * prime <= n);

	low3 = low_value; high3 = high_value; start3 = 0; count = 0;
	
	if (size % cacheSize == 0) end3 = size / cacheSize;
	
	else end3 = (size / cacheSize) + 1;

	do {
		low_value = low3 + (start3 * cacheSize * 2);
		high_value = MIN(high3, (low_value + (2 * cacheSize - 2)));
		size3 = (high_value - low_value) / 2 + 1;

		for (i = 0; i < size3; i++) marked[i] = 0;

		index = 0; prime = 3;

		do {
			if (currMarked[index] == 0) {
				prime = 2 * index + 3;
		
				if (prime * prime > low_value) {
					first = (prime * prime - low_value) / 2;
				}
				else {
					if (!(low_value % prime)) {
						first = 0;
					}
					else {
						if ((low_value % prime) % 2 == 0) {
							first = prime - (low_value % prime) / 2;
						} else {
							first = (prime - (low_value % prime)) / 2;
						}
					}
				}

				for (i = first; i < size3; i += prime) marked[i] = 1;
			}
			index++;
		} while (index < (((currHigh - 3) / 2) + 1));
		
		for (i = 0; i < size3; i++) {
			if (!marked[i]) count++;
		}
		start3++;
	} while (start3 < end3);

   	if (p > 1) {
		MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   	}

	elapsed_time += MPI_Wtime();
   	
	if (!id) {
		global_count++;
      		printf ("%llu primes are less than or equal to %llu\n", global_count, n);
      		printf ("Total elapsed time: %10.6f\n", elapsed_time);
   	}
   	
	MPI_Finalize();
   	
	return 0;
}
