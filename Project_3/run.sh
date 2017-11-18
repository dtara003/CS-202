#!/bin/sh

cd $PBS_O_WORKDIR

echo "Project 3: Parallel Sieve of Eratosthenes for Finding All Prime Numbers within 10^10"

module purge
module load gcc-4.6.2
module load mvapich2-1.9a2/gnu-4.6.2

mpicc -o part1 part1.c -lm
#mpicc -o part2 part2.c -lm
#mpicc -o part3 part3.c -lm

joba=$(qsub p11job)
#jobb=$(qsub -W depend=afterany:$(joba) p12job)
#jobc=$(qsub -W depend=afterany:$(jobb) p14job)
#jobd=$(qsub -W depend=afterany:$(jobc) p18job)
