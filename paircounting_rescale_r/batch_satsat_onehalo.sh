#!/bin/bash -l

#!/bin/bash
#SBATCH -p regular
#SBATCH -o logs/paircounting_test_satsat_onehalo
#SBATCH --time=40
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
# load nbodykit
#SBATCH --mail-type=END    # notifications for job
#SBATCH --mail-user=cameron.grove@durham.ac.uk

python3 satsat_onehalo.py


