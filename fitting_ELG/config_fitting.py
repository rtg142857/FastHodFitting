# config_fitting.py
"""
A file to contain input parameters for fitting functions
"""

import numpy as np
from scipy.special import erf

# Import values from paircounting config file
import config
import numpy as np
path = config.path

boxsize = config.boxsize

r_bin_edges = config.r_bin_edges

mass_bin_edges = config.mass_bin_edges

num_sat_parts = config.num_sat_parts

run_label = config.run_label


# Number of steps to take with each emcee walker
num_steps = 2000

# Number of walkers to use in the emcee fitting 
num_walkers = 20

# Run label for paircounting so the results can be read in
# Could make this dynamic so it used run label from above
run_path = "/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/ELG/high_res/FastHodFitting/paircounting/ELG_high"

# Include cen and sat HOD definitions in here as they can change when fitting different things

def Cen_HOD(p0,x):
    fa,fb,siga,sigb,Mc = p0[:5]
    x = np.log10(x)
    gau = np.zeros(shape=(len(x)))
    step = np.zeros(shape=(len(x)))

    ind = np.where(x<Mc)
    gau[ind] = fb*(1.-fa)/np.exp((x[ind]-Mc)**2./(2.*siga**2.))
    step[ind] = fa*(1.+erf((x[ind]-Mc)/siga))

    ind = np.where(x>=Mc)
    gau[ind] = fb*(1.-fa)/np.exp((x[ind]-Mc)**2./(2.*sigb**2.))
    step[ind] = fa*(1.+erf((x[ind]-Mc)/sigb))

    y = gau + step

    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    return 10**y


def Sat_HOD(p0,x):
    fs,Mmin,sigs,asat = p0[5:]
    x = np.log10(x)
    y = np.zeros(shape=(len(x)))
    if (Mmin>0. and sigs>0.):
        ind = np.where(x>0.) ; negs = np.shape(ind)[1]-np.shape(x)[0]
        if (negs>0):
            print('lSat_C13, WARNING: Negative masses as input')
        y[ind] = fs*(1.+erf((x[ind]-Mmin)/sigs))*\
            (10**(x[ind])/10**(Mmin))**asat

    ind = np.where(y>0.)
    if (np.shape(ind)[1]>1):
        y[ind] = np.log10(y[ind])
    y[x<10.8] = -100.
    return 10**y

# Likelihood definition can change a lot between fits as well 

def likelihood_calc(model,y,err):
    # Here take a constant fractional error and exclude BAO scale
    likelihood = - 0.5 * (np.sum((1 - (model[:] / y[:]))**2) / err**2)
    return likelihood
    
# Target correlation function to fit to
# Rescale this by the cosmology factor here:
#cosmo_factor = np.genfromtxt("cosmology_rescaling_factor.txt")
target_2pcf = np.genfromtxt("xir_Abacus-HOD-high_z1.dat")
#for i in range(10):
#    target_2pcf[:,i] = target_2pcf[:,i] * cosmo_factor

# Target number density array
target_num_den = 1e-4 #np.genfromtxt("/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/target_number_density.dat")
# error to apply to each point on the correlation function, can affect the speed of the fitting and walkers
# can get stuck in local minima if this is set too small (<0.1) 
err = .1

# Number density error parameter, affects how tight the target number density constraint is
num_den_err = 1e-3

# Prior_locations, positions of the priors on the parameters for fitting

priors = np.array([[-1,100],
                   [-1,100],
                   [-1,100],
                   [-1,100],
                   [-1,100],                  
                   [-1,100],
                   [-1,100],
                   [-1,100],
                   [-1,100]

])

# Flag to change whether initial parameters are in random locations across the prior space or not
initial_params_random = False

# If initial_params_random = 0 then use initial parameters close to these values:
initial_params = np.array([4.358e-03, 0.065, 0.116, 0.273, 11.533,2.331e-04, 11.260, 0.188, 0.701])


num_mass_bins_big = 20000

mass_function = np.genfromtxt("/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/mass_function.dat")

random_seed = 1234

save_path = "original_params"
