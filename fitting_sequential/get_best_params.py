import numpy as np
import sys
import emcee

import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.special import erf,erfc
import scipy.integrate as integrate
import time
import corner
import pickle

samplers_path = sys.argv[1]

def log_probability():
    """unpickling gets confused unless this is included, maybe just save the emcee outputs as param and likelihood arrays in the future"""
    return 0


print(samplers_path)

with open(samplers_path, 'rb') as input:
        samplers = pickle.load(input, encoding="latin1")

def max_like_params(sampler):
    """Get the parameters which provide the maximum likelihood from the sampling"""

    # Only take after the first 1000 steps to allow burn in
    flat_samples = sampler.backend.get_chain()[1000:,:,:]
    likelihoods = sampler.backend.get_log_prob()[1000:,:]

    # Find the maximum likelihood parameter position

    print(np.shape(likelihoods))
    best_param_index1 = np.argmax(likelihoods,axis = 0)
    best_param_index2 = np.argmax(sampler.backend.get_log_prob()[best_param_index1,np.arange(20)])
    best_params = flat_samples[best_param_index1[best_param_index2],best_param_index2,:]
    print(best_params)
    return best_params




best_fit_params = np.zeros((9,6))

number_density = np.genfromtxt("/cosma7/data/durham/dc-grov1/Halo_mass_pair_binning/BGS/target_number_density.dat")

best_fit_params[:,0] = number_density[:,0]

# Number of fitted params changes depending on how far through the parameter fixing we are, include an iteration number to solve this

iteration_number = int(sys.argv[2])

if iteration_number == 0:
    for i in range(9):
        best_fit_params[i,1:] = max_like_params(samplers[i])

elif iteration_number == 1:
    sigma_logM = np.genfromtxt("sigma_logM_fit.txt")
    for i in range(9):
        best_fit_params[i,[1,3,4,5]] = max_like_params(samplers[i])
    best_fit_params[:,2] = sigma_logM

elif iteration_number == 2:
    sigma_logM = np.genfromtxt("sigma_logM_fit.txt")
    alpha = np.genfromtxt("alpha_fit.txt")
    for i in range(9):
        best_fit_params[i,[1,3,4]] = max_like_params(samplers[i])
    best_fit_params[:,2] = sigma_logM
    best_fit_params[:,5] = alpha


elif iteration_number == 3:
    sigma_logM = np.genfromtxt("sigma_logM_fit.txt")
    alpha = np.genfromtxt("alpha_fit.txt")
    M0 = np.genfromtxt("M0_fit.txt")
    for i in range(9):
        best_fit_params[i,[1,4]] = max_like_params(samplers[i])
    best_fit_params[:,2] = sigma_logM
    best_fit_params[:,5] = alpha
    best_fit_params[:,3] = M0


elif iteration_number == 4:
    sigma_logM = np.genfromtxt("sigma_logM_fit.txt")
    alpha = np.genfromtxt("alpha_fit.txt")
    M0 = np.genfromtxt("M0_fit.txt")
    M1 = np.genfromtxt("M1_fit.txt")
    for i in range(9):
        best_fit_params[i,1] = max_like_params(samplers[i])
    best_fit_params[:,2] = sigma_logM
    best_fit_params[:,5] = alpha
    best_fit_params[:,3] = M0
    best_fit_params[:,4] = M1


np.savetxt("hod_params/params/fit_params.txt",best_fit_params)
















    
