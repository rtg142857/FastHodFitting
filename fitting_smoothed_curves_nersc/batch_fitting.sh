#!/bin/bash -l

#!/bin/bash
#SBATCH -p regular
#SBATCH -o logs/fitting_1
#SBATCH --time=600
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
# load nbodykit
#SBATCH --mail-type=END    # notifications for job
#SBATCH --mail-user=cameron.grove@durham.ac.uk




python3 fasthod_fitting.py 1



