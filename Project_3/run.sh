#!/bin/sh

cd $PBS_O_WORKDIR

echo "Project 3: Parallel Sieve of Eratosthenes for Finding All Prime Numbers within 10^10" 

module purge
module load gcc-4.6.2
module load mvapich2-1.9a2/gnu-4.6.2

mpicc -o p1 p1mpi.c -lm
mpicc -o p2 p2mpi.c -lm
mpicc -o p3 p3mpi.c -lm
  
jobid=$(qsub jobfile_part1_node1)
joba=$(qsub -W depend=afterany:${jobid} jobfile_part1_node2)
jobb=$(qsub -W depend=afterany:${joba} jobfile_part1_node4)
jobc=$(qsub -W depend=afterany:${jobb} jobfile_part1_node8)

jobd=$(qsub -W depend=afterany:${jobc} jobfile_part2_node1)
jobe=$(qsub -W depend=afterany:${jobd} jobfile_part2_node2)
jobf=$(qsub -W depend=afterany:${jobe} jobfile_part2_node4)
jobg=$(qsub -W depend=afterany:${jobf} jobfile_part2_node8)

jobh=$(qsub -W depend=afterany:${jobg} jobfile_part3_node1)
jobi=$(qsub -W depend=afterany:${jobj} jobfile_part3_node2)
jobj=$(qsub -W depend=afterany:${jobi} jobfile_part3_node4)
qsub -W depend=afterany:${jobj} jobfile_part3_node8
