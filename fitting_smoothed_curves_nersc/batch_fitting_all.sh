#!/bin/bash -l

#!/bin/bash
#SBATCH -p regular
#SBATCH -o logs/fitting_1
#SBATCH --time=1800
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
# load nbodykit
#SBATCH --mail-type=END    # notifications for job
#SBATCH --mail-user=cameron.grove@durham.ac.uk



python3 fasthod_fitting.py 1
python3 fasthod_fitting.py 2
python3 fasthod_fitting.py 3
python3 fasthod_fitting.py 4
python3 fasthod_fitting.py 5
python3 fasthod_fitting.py 6
python3 fasthod_fitting.py 7
python3 fasthod_fitting.py 8
python3 fasthod_fitting.py 9
python3 fasthod_fitting.py 10


