#!/bin/bash -l

#!/bin/bash
#SBATCH -q regular
#SBATCH -o logs/tracer_snap
#SBATCH --time=720
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH -C cpu

# activate DESI environment on NERSC
cosmodesienv main

module load gcc
module load gsl
module unload craype-hugepages2M


# Only do the halo paircounting in this step
python cencen.py
python censat.py
python satsat.py
python satsat_onehalo.py

