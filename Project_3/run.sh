#!/bin/sh
PBS_O_WORKDIR="/home/dtara003/CS-202/Project_3"

cd $PBS_O_WORKDIR

echo "Project 3: Parallel Sieve of Eratosthenes for Finding All Prime Numbers within 10^10"

module purge
module load gcc-4.6.2
module load mvapich2-1.9a2/gnu-4.6.2

mpicc -o part1 part1.c -lm
#mpicc -o part2 part2.c -lm
#mpicc -o part3 part3.c -lm

jobid=$(qsub p11job)
joba=$(qsub -W depend=afterany:${jobid} p12job)
jobb=$(qsub -W depend=afterany:${joba} p14job)
jobc=$(qsub -W depend=afterany:${jobb} p18job)

