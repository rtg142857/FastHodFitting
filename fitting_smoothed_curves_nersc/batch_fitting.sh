#!/bin/bash -l

#!/bin/bash
#SBATCH -q regular
#SBATCH -o log_fitting
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH -C cpu
#SBATCH --mail-user=alexander.m.smith@durham.ac.uk 

cosmodesienv main

module load gcc
module load gsl
module unload craype-hugepages2M

python3 fasthod_fitting.py 1
