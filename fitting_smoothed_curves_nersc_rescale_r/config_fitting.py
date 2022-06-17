# config_fitting.py
"""
A file to contain input parameters for fitting functions
"""

import numpy as np
from scipy.special import erf

# Import values from paircounting config file
import config
import numpy as np


rescale_r = 0.9758508296342128



path = config.path

boxsize = config.boxsize

r_bin_edges = config.r_bin_edges

mass_bin_edges = config.mass_bin_edges

num_sat_parts = config.num_sat_parts

run_label = config.run_label


# Number of steps to take with each emcee walker
num_steps = 2000

# Number of walkers to use in the emcee fitting 
num_walkers = 180

# Run label for paircounting so the results can be read in
# Could make this dynamic so it used run label from above
run_path = "../paircounting_rescale_r/pairs"

# Include cen and sat HOD definitions in here as they can change when fitting different things
def spline_kernel_integral(x):
    """
    Returns the integral of the unscaled spline kernel function from -1 to x
    """
    if hasattr(x, "__len__"):
        # x in an array
        integral = np.zeros(len(x))
        absx = abs(x)
        ind = absx < 0.5
        integral[ind] = absx[ind] - 2*absx[ind]**3 + 1.5*absx[ind]**4
        ind = np.logical_and(absx >= 0.5, absx < 1.)
        integral[ind] = 0.375 - 0.5*(1-absx[ind])**4
        ind = absx >= 1.
        integral[ind] = 0.375
        ind = x < 0
        integral[ind] = -integral[ind]
    else:
        # x is a number
        absx = abs(x)
        if   absx < 0.5: integral = absx - 2*absx**3 + 1.5*absx**4
        elif absx < 1:   integral = 0.375 - 0.5*(1-absx)**4
        else:            integral = 0.375
        if x < 0: integral = -integral
    return integral

def cumulative_spline_kernel(x, mean=0, sig=1):
    """
    Returns the integral of the rescaled spline kernel function from -inf to x.
    The spline kernel is rescaled to have the specified mean and standard
    deviation, and is normalized.
    """
    integral = spline_kernel_integral((x-mean)/(sig*np.sqrt(12))) / 0.75
    y = 0.5 * (1. + 2*integral)
    return y

# Now change Cen_HOD definition

def Cen_HOD(params,mass_bins):
    Mmin, sigma_logm = params[:2]
    result = cumulative_spline_kernel(np.log10(mass_bins), mean = Mmin, sig=sigma_logm/np.sqrt(2))
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
    #likelihood = - 0.5 * (np.sum((1 - (model[:75] / y[:75]))**2) / err**2)
<<<<<<< HEAD
    # Alternatively use a less tight fit at very small scales, using observed effect of sample variance
    # cutoff at large scales so systematic effects aren't overweighted
=======
    # A temporary change to test a feature:
>>>>>>> 256e5bd3e652318599aa8a60c089b4509bc16e08
    likelihood = - 0.5 * (np.sum((1 - (model[:75] / y[:75]))**2 / (2*np.maximum(0.01 * np.logspace(-2,2.2,85)[:75]**-0.5,0.03))**2))
    return likelihood
    
# Target correlation function to fit to
# Rescale this by the cosmology factor here:
cosmo_factor = np.genfromtxt("cosmology_rescaling_factor_zel_8.txt")
target_2pcf = np.genfromtxt("xi_r_mxxl.dat")
for i in range(10):
    target_2pcf[:,i] = target_2pcf[:,i] * cosmo_factor

# Target number density array
<<<<<<< HEAD
target_num_den = np.genfromtxt("target_num_den_rescaled.dat")
=======
# target_num_den = np.genfromtxt("target_number_density.dat")
target_num_den = np.genfromtxt("target_num_den_c004.txt") # Temporary ajustment for test
>>>>>>> 256e5bd3e652318599aa8a60c089b4509bc16e08

# error to apply to each point on the correlation function, can affect the speed of the fitting and walkers
# can get stuck in local minima if this is set too small (<0.1) 
err = 0.7 * 0.0625 * 1.5

# Number density error parameter, affects how tight the target number density constraint is
num_den_err = 0.0625 * (0.01**0.5) * 1.5

# Prior_locations, positions of the priors on the parameters for fitting

priors = np.array([[10,16],
                   [0,4],
                   [0,16], # Changed
                   [10,17],
                   [0,4]
])

# Add in meta-priors for the parameters defining the shape of the smooth curves

meta_priors = np.array([[-5,5],
                   [-3,3],
                   [-3,3], # Mmin param 3
                   [-3,3],
                   [2e-3,2e-1],
                   [0.06,6.0],
                   [18.0,24.0],
                   [0.25,25.0],
                   [-10,10], # Changed
                   [-20,20], # Changed
                   [-5.,5.],
                   [-3.,3.],
                   [-3.,3.], # M1 param3
                   [-3.,3.0],
                   [0.4,1.5],
                   [2.,15.0],
                   [-5.,0.0]
])


# Flag to change whether initial parameters are in random locations across the prior space or not
initial_params_random = False

# If initial_params_random = 0 then use initial parameters close to these values:
#initial_params = np.array([1, -0.1,0.01, -0.001,#Mmin
#0.03, 0.8, 20, 4,#sigmalogm
#-0.2, 0.4,#M0
#1, -0.1, 0.01,-0.001,#M1
#1, 9, -3])#alpha

initial_params = np.array([
-1.103954732697280811e-01,
-5.607159196804438750e-01,
1.335559408347629373e-01,
-1.727281984610109372e-02,
1.547358691854127943e-01,
6.797987495699331362e-01,
2.025631753717012984e+01,
1.239389028687060357e+00,
-2.083123726627407635e+00,
-1.533159445569658885e+00,
1.083286411184483100e+00,
-4.839587028987296091e-01,
9.402281378285631819e-02,
-1.023851152682873390e-02,
1.095817887824390935e+00,
3.983897146373905684e+00,
-3.417082491894401830e+00
])

num_mass_bins_big = 18000


random_seed = 10

save_path = "fits_test_wider_start1"
