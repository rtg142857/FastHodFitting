#!/bin/bash -l

#!/bin/bash
#SBATCH -p regular
#SBATCH -o logs/paircounting_test_censat
#SBATCH --time=100
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
# load nbodykit
#SBATCH --mail-type=END    # notifications for job
#SBATCH --mail-user=cameron.grove@durham.ac.uk

module load gcc
module load gsl
module unload craype-hugepages2M
python3 censat.py


