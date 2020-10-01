# config_fitting.py
"""
A file to contain input parameters for fitting functions
"""

import numpy as np
from scipy.special import erf

# Import values from paircounting config file
import config

path = config.path

boxsize = config.boxsize

r_bin_edges = config.r_bin_edges

mass_bin_edges = config.mass_bin_edges

num_sat_parts = config.num_sat_parts

run_label = config.run_label


# Number of steps to take with each emcee walker
num_steps = 1000

# Number of walkers to use in the emcee fitting 
num_walkers = 20

# Run label for paircounting so the results can be read in
# Could make this dynamic so it used run label from above
run_path = "/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/Github/FastHodFitting/paircounting/test"

# Include cen and sat HOD definitions in here as they can change when fitting different things

def Cen_HOD(params,mass_bins):
    Mmin, sigma_logm = params[:2].copy()
    Mmin = 10**Mmin
    result = 0.5*(1+erf((np.log10(mass_bins)-np.log10(Mmin))/sigma_logm))
    return(result)

def Sat_HOD(params,cen_hod,mass_bins):
    M0, M1, alpha = params[2:].copy()
    M0 = 10**M0
    M1 = 10**M1
    result = cen_hod * (((mass_bins-M0)/M1)**alpha)
    return(result)

# Likelihood definition can change a lot between fits as well 

def likelihood_calc(model,y,err):
    # Here take a constant fractional error and exclude BAO scale
    likelihood = - 0.5 * (np.sum((1 - (model[:75] / y[:75]))**2) / err**2)
    return likelihood
    
# Target correlation function to fit to
target_2pcf = np.genfromtxt("/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/xi_r_mxxl.dat")

# Target number density array
target_num_den = np.genfromtxt("/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/target_number_density.dat")

# error to apply to each point on the correlation function, can affect the speed of the fitting and walkers
# can get stuck in local minima if this is set too small (<0.1) 
err = 1.

# Number density error parameter, affects how tight the target number density constraint is
num_den_err = 0.01

# Prior_locations, positions of the priors on the parameters for fitting

priors = np.array([[11,16],
                   [0,4],
                   [8,14],
                   [11,17],
                   [0,4]
])

# Flag to change whether initial parameters are in random locations across the prior space or not
initial_params_random = False

# If initial_params_random = 0 then use initial parameters close to these values:
initial_params = np.array([12,0.1,10,13,0.9])

num_mass_bins_big = 20000

mass_function = np.genfromtxt("/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/mass_function.dat")

random_seed = 1234

save_path = "fits"
