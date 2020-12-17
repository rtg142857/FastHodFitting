#!/bin/bash -l

#SBATCH --ntasks 1
#SBATCH -J abacus_fits_corrected_rescaled_pg
#SBATCH -o ./logs/%x.%J.out
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH -t 12:00:00
#SBATCH --mail-type=END    # notifications for job
#SBATCH --mail-user=cameron.grove@durham.ac.uk


module purge
module load gnu_comp openmpi python/3.6.5

sed -i '31s/.*/included_params = [1,1,1,1,1]/' config_fitting.py

python3 fasthod_fitting.py
cd hod_params
python3 get_fits.py
cd ..
python3 get_fitted_param.py 2 # Fix sigma_logM

# Now do fits with no sigma_logM
sed -i '31s/.*/included_params = [1,0,1,1,1]/' config_fitting.py

python3 fasthod_fitting.py
cd hod_params
python3 get_fits.py
cd ..
python3 get_fitted_param.py 4 # Fix M1

sed -i '31s/.*/included_params = [1,0,1,0,1]/' config_fitting.py

python3 fasthod_fitting.py
cd hod_params
python3 get_fits.py
cd ..
python3 get_fitted_param.py 5 # Fix alpha

sed -i '31s/.*/included_params = [1,0,1,0,0]/' config_fitting.py

python3 fasthod_fitting.py
cd hod_params
python3 get_fits.py
cd ..
python3 get_fitted_param.py 3 # Fix M0

sed -i '31s/.*/included_params = [1,0,0,0,0]/' config_fitting.py

python3 fasthod_fitting.py
cd hod_params
python3 get_fits.py
cd ..
python3 get_fitted_param.py 1 #  Fix Mmin (final parameter)

# Finally rescale Mmin to get the exactly correct number density using Alex's mass function
# The final parameter values will be produced at: number_density/params_fixed_number_density.txt
cd number_density

python3 number_density_all_mags.py


