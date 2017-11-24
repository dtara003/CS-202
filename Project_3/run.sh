#!/bin/sh
PBS_O_WORKDIR="/home/dtara003/CS-202/Project_3"

cd $PBS_O_WORKDIR

echo "Project 3: Parallel Sieve of Eratosthenes for Finding All Prime Numbers within 10^10"

module purge
module load gcc-4.6.2
module load mvapich2-1.9a2/gnu-4.6.2

mpicc -o part1 part1.c -lm
mpicc -o part2 part2.c -lm
mpicc -o part3 part3.c -lm

jobid1=$(qsub p11job)
joba1=$(qsub -W depend=afterany:${jobid1} p12job)
jobb1=$(qsub -W depend=afterany:${joba1} p14job)
jobc1=$(qsub -W depend=afterany:${jobb1} p18job)

jobid2=$(qsub -W depend=afterany:${jobc1} p21job)
joba2=$(qsub -W depend=afterany:${jobid2} p22job)
jobb2=$(qsub -W depend=afterany:${joba2} p24job)
jobc2=$(qsub -W depend=afterany:${jobb2} p28job)

jobid3=$(qsub -W depend=afterany:${jobc2} p31job)
joba3=$(qsub -W depend=afterany:${jobid3} p32job)
jobb3=$(qsub -W depend=afterany:${joba3} p34job)
jobc3=$(qsub -W depend=afterany:${jobb3} p38job)

