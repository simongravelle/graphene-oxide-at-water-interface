#!/bin/sh 
#$ -cwd 
#$ -pe smp 4
#$ -l h_rt=100:0:0 
#$ -l h_vmem=1G 

module load gcc
module load openmpi

mpirun -np 4 /data/home/exx356/lammps-29Oct20/src/lmp_mpi -in input.lammps
