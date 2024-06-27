#!/bin/bash -l

# Command line argument is path to path_config.yml

#!/bin/bash
#SBATCH -J HOD_fitting
#SBATCH -o logs/fitting
#SBATCH -e logs/fitting_error
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH -p cosma8 #or some other partition, e.g. cosma, cosma8, etc.
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH --mail-user=tlrt88@durham.ac.uk 

module purge
module use /cosma/home/dp004/dc-mene1/software/desi/cosmodesiconda/my-desiconda/modulefiles
module load cosmodesiconda/my-desiconda
#cosmodesienv main

#module load gcc
module load gsl
#module unload craype-hugepages2M
module load python/3.10.12

python fasthod_fitting.py $1
