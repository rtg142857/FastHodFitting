# fasthod_fitting.py
"""
A module for all the functions used in my code to fit HOD parameters
This should be imported at the start of any scripts
"""

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import h5py
import sys
import time
import emcee
from multiprocessing import Pool
import config
import config_fitting
# Load parameters from config file
path = config.path
path_tracer = config.path_tracer

boxsize = config.boxsize
r_bin_edges = config.r_bin_edges
mass_bin_edges = config.mass_bin_edges
num_sat_parts = config.num_sat_parts
run_label = config.run_label
subsample_array = config.subsample_array
# Load parameters from fitting config file

num_steps = config_fitting.num_steps
num_walkers = config_fitting.num_walkers
run_path = config_fitting.run_path
Cen_HOD = config_fitting.Cen_HOD
spline_kernel_integral = config_fitting.spline_kernel_integral
cumulative_spline_kernel = config_fitting.cumulative_spline_kernel
Sat_HOD = config_fitting.Sat_HOD
likelihood_calc = config_fitting.likelihood_calc
target_2pcf = config_fitting.target_2pcf
target_num_den = config_fitting.target_num_den
err = config_fitting.err
num_den_err = config_fitting.num_den_err
priors = config_fitting.priors
initial_params_random = config_fitting.initial_params_random
initial_params = config_fitting.initial_params
num_mass_bins_big = config_fitting.num_mass_bins_big
random_seed = config_fitting.random_seed

save_path = config_fitting.save_path

# Load smoothed HOD parameter curves from file and add meta-priors

meta_priors = config_fitting.meta_priors


random_seed = int(sys.argv[1])
np.random.seed(random_seed)

initial_params = initial_params * np.random.uniform(size = 17,low = 0.9, high = 1.1)
save_path = "diff_start_low_prior_"+str(random_seed)
print(initial_params)
print(save_path)
# quickly adjust initial params to get a random new starting position


def M0_function(magnitude, A, B):
    M0s = (A*(magnitude+20) + (B+11))
    M0s[M0s <=1.0] = 1.0
    return M0s

def sigma_function(magnitude, A, B, C, D):
    return A + (B-A) / (1.+np.exp((magnitude+C)*D))


def alpha_function(magnitude,A,B,C):
    return A + B ** (-magnitude - 20 + C)


def L_function(magnitude,A,B,C,D):
    return (A + 12) + B*(magnitude + 20) + C*(magnitude+20)**2 + D*(magnitude+20)**3 



def initialise_walkers(initial_params_random,initial_params,priors,num_walkers):
    """
    Initialise the positions of the walkers for fitting the HOD parameters
    Do this randomly within the prior space if initial_params_random=True
    Else populate in a small region around some provided params
    """
    pos = np.zeros((num_walkers,np.shape(initial_params)[0]))
    if (initial_params_random):
        for i in range(num_walkers):
            for j in range(np.shape(meta_priors)[0]):
                pos[i,j] = np.random.uniform(meta_priors[j,0],meta_priors[j,1])

    else:
        #if len(initial_params)!=np.shape(priors)[0]:
        #    raise ValueError("Your initial parameter values and priors have different shapes")
        for i in range(num_walkers):
            # Populate in a 10% region around parameters provided
            # Potential to make the size of this region an input variable if necessary
            pos[i,:] = initial_params*(0.95 + 0.1*np.random.random(np.shape(initial_params)[0]))
    #print(pos)

    # Check none of the walkers lie outside the prior space
    #for i in range(num_walkers):
    #    for j in range(np.shape(priors)[0]):
    #        if pos[i,j] < priors[j,0]:
    #            raise ValueError("Your initial parameter values lie outside the prior space, parameter ",j, " is too low")
    #        if pos[i,j] > priors[j,1]:
    #            raise ValueError("Your initial parameter values lie outside the prior space, parameter ",j, " is too high")
    print(pos)
    return pos

        
def undo_subsampling(mass_pair_array,subsample_array):
    """
    Change a halo paircount table to undo the effect of subsampling
    """
    altered_mass_pair_array = np.zeros(np.shape(mass_pair_array))
    for i in range(len(mass_pair_array[:,0,0])):
        for j in range(len(mass_pair_array[:,0,0])):
            altered_mass_pair_array[i,j,:] = mass_pair_array[i,j,:] / (subsample_array[i]*subsample_array[j])
    return altered_mass_pair_array 

def calc_hmf(path,num_mass_bins_big,mass_bin_edges):
    """
    Calculate the hmf from the catalog
    """
    snap = h5py.File(path,"r")
    Mvir = snap["/mass"][:]*1e10
    is_central = snap["/is_central"][:]


    mass_centrals = Mvir[is_central]
    mass_sats = Mvir[~np.array(is_central)]

    tracer_snap = h5py.File(path_tracer,"r")

    mass_centrals = np.append(mass_centrals,tracer_snap["/mass"][:]*1e10)

    mass_min = mass_bin_edges[0]
    mass_max = mass_bin_edges[-1]
    
    
    # Go for a large number of sub bins for accuracy
    mass_bins_big = np.logspace(np.log10(mass_min),np.log10(mass_max),num_mass_bins_big + 1)
    cen_halos_big = np.histogram(mass_centrals,bins = mass_bins_big)[0]
    sat_halos_big = np.histogram(mass_sats,bins = mass_bins_big)[0]
    
    return mass_bins_big, cen_halos_big, sat_halos_big

def calc_hmf_more_files(path,num_mass_bins_big,mass_bin_edges):
    """
    Calculate the hmf from the catalog
    New file format split into 34 separately
    """
    snap = h5py.File(path+"galaxy_tracers_0.hdf5","r")
    Mvir = snap["/mass"][:]
    is_central = snap["/is_central"][:]

    # There won't be double the number of particles
    # in another file compared to the first one
    # So use this to create unique halo IDs
    for i in range(1,34):
        snap = h5py.File(path+"galaxy_tracers_"+str(i)+".hdf5","r")
        Mvir_temp = snap["/mass"][:]
        is_central_temp = snap["/is_central"][:]

        Mvir = np.append(Mvir,Mvir_temp)
        is_central = np.append(is_central,is_central_temp)
        print("Reading File Number "+str(i))


    mass_centrals = Mvir[is_central]*1e10
    mass_sats = Mvir[~np.array(is_central)]*1e10

    snap = h5py.File(path+"galaxy_tracers_unresolved_0.hdf5","r")
    Mvir = snap["/mass"][:]

    # There won't be double the number of particles
    # in another file compared to the first one
    # So use this to create unique halo IDs
    for i in range(1,34):
        snap = h5py.File(path+"galaxy_tracers_unresolved_"+str(i)+".hdf5","r")
        Mvir_temp = snap["/mass"][:]

        Mvir = np.append(Mvir,Mvir_temp)
        print("Reading Unresolved Halo File Number "+str(i))
    mass_centrals = np.append(mass_centrals,Mvir*1e10)

    mass_min = mass_bin_edges[0]
    mass_max = mass_bin_edges[-1]


    # Go for a large number of sub bins for accuracy
    mass_bins_big = np.logspace(np.log10(mass_min),np.log10(mass_max),num_mass_bins_big + 1)
    cen_halos_big = np.histogram(mass_centrals,bins = mass_bins_big)[0]
    sat_halos_big = np.histogram(mass_sats,bins = mass_bins_big)[0]

    return mass_bins_big, cen_halos_big, sat_halos_big

def calc_hmf_ascii(path, num_mass_bins_big,mass_bin_edges):
    """
    Calculate the hmf from the catalog
    """

    path = sys.argv[1]
    Mvir, PID = np.loadtxt(path,usecols=(2,41),unpack=True)


    # Only take the central halos as we shall create our own satellite particles later
    mask1 = np.where(PID == -1)

    Mvir = Mvir[mask1]


    mass_min = mass_bin_edges[0]
    mass_max = mass_bin_edges[-1]


    # Go for a large number of sub bins for accuracy
    mass_bins_big = np.logspace(np.log10(mass_min),np.log10(mass_max),num_mass_bins_big + 1)
    cen_halos_big = np.histogram(Mvir,bins = mass_bins_big)[0]
    sat_halos_big = np.histogram(Mvir,bins = mass_bins_big)[0]

    return mass_bins_big, cen_halos_big, sat_halos_big
 


  
def create_weighting_factor(mass_pair_array,hod1,hod2):
    """
    Multiply the array by the relevant HODs and then sum over the mass bins
    to get the number of pairs as a function of r. These can then be divided
    by the randoms to get the correlation function.
    """
    weighting_factor = np.tensordot(np.outer(hod1,hod2),mass_pair_array,axes=([0,1],[0,1]))
    return weighting_factor

def create_accurate_HOD(hod,halos,mass_bin_edges,num_mass_bins_big):
    """
    Using just 100-400 mass bins for the HOD isn't accurate enough. Take much smaller
    mass subdivisions and use these to create an accurate HOD for only 100-400 mass bins.
    """
    num_mass_bins_small = int(len(mass_bin_edges) -1)
    mass_bins_factor = int(num_mass_bins_big/num_mass_bins_small)
    if num_mass_bins_big % num_mass_bins_small != 0:
        raise ValueError("finer grained mass bins do not evenly divide coarser mass bins:",num_mass_bins_small," is not a factor of ",num_mass_bins_big)
    hod[np.isnan(hod)] = 0
    HOD_halo_product = halos * hod
    HOD_recalc = np.sum(np.reshape(HOD_halo_product,(num_mass_bins_small,mass_bins_factor)),axis=1) / (
                 np.sum(np.reshape(halos,(num_mass_bins_small,mass_bins_factor)),axis=1))
    HOD_recalc[np.isnan(HOD_recalc)] = 0
    return HOD_recalc


def create_2pcf(hod_cen,hod_sat,cencen,censat,satsat,satsat_onehalo,num_sat_parts,
                hod_cen_big,hod_sat_big,cen_halos_big,sat_halos_big,
                boxsize):
    """
    Creates the 2pcf from the hod and the paircounts.
    First create the total number of paircounts.
    Then create the analytic randoms based on the total galaxies expected
    Then create the correlation function from these two components
    """
    
    weighting_factor_cencen = create_weighting_factor(cencen,hod_cen,hod_cen)
    # Multiply censat factor by two as the others are double counted but this one isn't
    weighting_factor_censat = create_weighting_factor(censat,hod_cen,hod_sat)*2
    weighting_factor_satsat = create_weighting_factor(satsat,hod_sat,hod_sat)
    weighting_factor_satsat_onehalo = create_weighting_factor(satsat_onehalo,hod_sat,hod_sat) / (
                                      (num_sat_parts*(num_sat_parts-1))/2)

    weighting_factor_total = weighting_factor_cencen + weighting_factor_censat + (
                             weighting_factor_satsat + weighting_factor_satsat_onehalo)

    npart = np.sum(cen_halos_big * hod_cen_big)
    npart_sat = np.sum(sat_halos_big * hod_sat_big / num_sat_parts)
    randoms = create_randoms_for_corrfunc(boxsize = boxsize,
                                          npart = (npart + npart_sat),
                                          r_bin_edges = r_bin_edges)
    cf = (weighting_factor_total / randoms) - 1
    return cf



def create_randoms_for_corrfunc(npart,r_bin_edges,boxsize):
    """
    Calculate the analytic randoms for npart particles in a box with
    side length boxsize.  This code is based on the calculation done 
    either in halotools or corrfunc but the formula is pretty simple.
    """
    NR = npart

    # do volume calculations
    v = (4./3.)*np.pi*r_bin_edges**3 # volume of spheres
    dv = np.diff(v)  # volume of shells
    global_volume = boxsize**3  # volume of simulation

    # calculate the random-random pairs using density * shell volume
    rhor = (NR*(NR-1))/global_volume
    RR = (dv*rhor)
    return RR

def calc_number_density(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize):
    """
    A function which caluculates the number density of halos from the hod
    and a mass function
    """
    num_density = (np.sum(hod_cen_big*cen_halos_big) + np.sum(hod_sat_big*cen_halos_big)) / (boxsize**3)

    return num_density
    
def create_hod_params(params):
    """
    Create the HOD parameters from the smoothed curve meta-parameters
    Smoothed curve functions can be defined at the top of this script.
    Returns arrays containing HOD parameters at each magnitude limit
    in the samples
    """
    mags = target_num_den[:,0]
    MminA,MminB,MminC,MminD, slmA, slmB, slmC, slmD, M0A, M0B, M1A, M1B, M1C,M1D, alA, alB, alC = params
    Mmin = L_function(mags,MminA,MminB,MminC,MminD)
    slm = sigma_function(mags,slmA,slmB,slmC,slmD)
    M0 = M0_function(mags,M0A,M0B)
    M1 = L_function(mags,M1A,M1B,M1C,M1D)
    alpha = alpha_function(mags,alA,alB,alC)
    HOD_params = np.zeros((len(Mmin),5))
    HOD_params[:,0] = Mmin
    HOD_params[:,1] = slm
    HOD_params[:,2] = M0
    HOD_params[:,3] = M1
    HOD_params[:,4] = alpha
    #print("Mmin",Mmin[0],"slm",slm[0],"M0",M0[0],"M1",M1[0],"alpha",alpha[0])
    return HOD_params

def log_likelihood_calc(params, target, target_density,target_number):
    """
    Calculates the contribution to the likelihood from an individual
    sample at one magnitude limit. The total likelihood function
    adds up the total contribution from all samples
    """
    hod_cen_big = Cen_HOD(params,mass_bin_centres_big)
    hod_sat_big = Sat_HOD(params,hod_cen_big,mass_bin_centres_big)
    
    
    hod_cen = create_accurate_HOD(hod_cen_big,cen_halos_big,mass_bin_edges,num_mass_bins_big)
    hod_sat = create_accurate_HOD(hod_sat_big,sat_halos_big,mass_bin_edges,num_mass_bins_big)
    cf = create_2pcf(hod_cen,hod_sat,cencen,censat,satsat,satsat_onehalo,num_sat_parts,
                hod_cen_big,hod_sat_big,cen_halos_big,sat_halos_big,
                boxsize)

    number_density = calc_number_density(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize)

    # Likelihood from fitting cf
    if target_number <=8:
        likelihood = likelihood_calc(cf,target[:,target_number+1],err)
    else:
        likelihood = 0

    # Likelihood from number density fit
    likelihood_num_den =  - (1 - (number_density / target_density))**2 / (num_den_err**2)

    # Add them both together to get the total likelihood
    l_likelihood = likelihood + likelihood_num_den
    return l_likelihood


def log_likelihood(params, target, target_density):
    """
    A function to convert from the parameters to the log-likelihood
    Involves creating the 2pcf, then comparing it to the model
    Also includes matching the number density as well
    """
    hod_params = create_hod_params(params)
    like_total = 0
    for i in range(len(hod_params[:,0])):
        like_total += log_likelihood_calc(hod_params[i,:],target,target_density[i],i)
    l_likelihood = like_total / (len(hod_params[:,0])-2) # To move to the same scale as previous runs, -2 is for the -17,-17.5 targets
    return l_likelihood
    
    
def log_prior_calc(theta):
    outside_prior_counter = 0
    for i in range(np.shape(priors)[0]):
        if ((theta[i] < priors[i,0]) or (theta[i] > priors[i,1]) ):
            outside_prior_counter +=1

    return outside_prior_counter

    

def log_prior(theta):
    """
    Returns whether the parameters are within the priors or not
    If parameters are outside any priors then return a -inf log probability
    """
    hod_params = create_hod_params(theta)
    outside_prior_counter = 0
    for i in range(len(hod_params[:,0])):
        outside_prior_counter += log_prior_calc(hod_params[i,:])

    for i in range(np.shape(meta_priors)[0]):
        if ((theta[i] < meta_priors[i,0]) or (theta[i] > meta_priors[i,1]) ):
            outside_prior_counter +=1

    if outside_prior_counter == 0:
        return 0.0
    else:
        return -np.inf


def log_probability(theta, y, yerr):
    """
    Converts the parameter values, target 2pcf, and error into a probability for emcee to use
    """
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(theta, y, yerr)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll


def perform_fitting(target_number):
    """
    A function to perform the fitting of the parameters
    Takes only target number as an input, this is the index
    of the magnitude limit which is being fitted. This number
    will be iterated over
    """
    walker_init_pos = initialise_walkers(initial_params_random,initial_params,priors,num_walkers)
    
    target = target_2pcf[:,:]
    target_density = 10**target_num_den[:,1]

    #target = target_2pcf[:,target_number+1]
    #target_density = 10**target_num_den[target_number,1]
    print("Luminosity limit: ",target_num_den[target_number,0])
    print("Target Density: ",10**target_num_den[target_number,1])
    nwalkers, ndim = walker_init_pos.shape
    start_time = time.time()
    #with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(target, target_density))
    sampler.run_mcmc(walker_init_pos, num_steps);
    end_time = time.time()
    print("fitting took ", end_time - start_time, " seconds")
    return sampler

def plot_HODs(HODs):
    """
    Plot HOD/s resulting from fits
    """
    plt.figure(figsize = (8,8))
    for i in range(len(HODs[0,:])):
        plt.plot(mass_bin_edges[:-1] + np.diff(mass_bin_edges)/2,HODs[:,i],c="C"+str(i),label=target_num_den[i,0])
    plt.ylim(1e-3)
    plt.yscale("log")
    plt.xscale("log")
    plt.legend(title = "Magnitude")
    plt.ylabel("n")
    plt.xlabel("Halo Mass /Solar Masses")
    plt.savefig(save_path+"_HODs.png",bbox_inches="tight")
    plt.close()
    return 0

def plot_CFs(CFs):
    """
    Plot CF ratio/s from fits
    """
    #plt.figure(figsize = (10,8))
    for i in range(len(CFs[0,:])):
        r_bin_centres = 10**(np.log10(r_bin_edges[:-1]) + np.diff(np.log10(r_bin_edges))/2 )
        plt.plot(r_bin_centres,CFs[:,i],c="C"+str(i))#,linestyle="--",alpha=0.6)
    plt.xscale("log")
    plt.xlabel("r [Mpc/h]")
    plt.ylabel(r"$\xi$(r) ratio")
    plt.legend(title = "Magnitude")
    plt.ylim(0.7,1.5)
    plt.grid()
    plt.title("Log Likelihood: "+str(np.round(best_log_prob,4)))
    plt.savefig(save_path+"_corrfunc_ratio.png",bbox_inches="tight")
    plt.close()
    return 0

def plot_CFs_errors(CFs):
    """
    Plot CF ratio/s from fits
    """
    #plt.figure(figsize = (10,8))
    for i in range(len(CFs[0,:])):
        r_bin_centres = 10**(np.log10(r_bin_edges[:-1]) + np.diff(np.log10(r_bin_edges))/2 )
        plt.plot(r_bin_centres,CFs[:,i],c="C"+str(i))#,linestyle="--",alpha=0.6)
    plt.hlines(xmin = r_bin_centres[0],xmax =r_bin_centres[74],y=1+ 3 * 0.7 * 0.0625 * 1.5,colors="k")
    plt.hlines(xmin = r_bin_centres[0],xmax =r_bin_centres[74],y=1- 3 * 0.7 * 0.0625 * 1.5,colors="k")
    plt.xscale("log")
    plt.xlabel("r [Mpc/h]")
    plt.ylabel(r"$\xi$(r) ratio")
    plt.legend(title = "Magnitude")
    plt.ylim(0.7,1.5)
    plt.grid()
    #plt.title("Log Likelihood: "+str(np.round(best_log_prob,4)))
    plt.savefig(save_path+"_corrfunc_ratio_errors.png",bbox_inches="tight")
    plt.close()
    return 0


def plot_HOD_params(HOD_params):
    """
    Plot HOD parameters against magnitude
    """
    
    labels = ["Mmin","sigma_logm","M0","M1","alpha"]
    fig, axs = plt.subplots(ncols=2,nrows=3,figsize=(7,10))
    for i in range(5):
        axs[int(np.floor(i/2)),i%2].errorbar(x=-target_num_den[:,0],y=HOD_params[:,i],yerr = 0,color="C0")
        axs[int(np.floor(i/2)),i%2].set_xlabel("Magnitude")
        axs[int(np.floor(i/2)),i%2].set_ylabel("Parameter Value")
        axs[int(np.floor(i/2)),i%2].legend(title=labels[i])
    i = 5
    axs[int(np.floor(i/2)),i%2].errorbar(x=-target_num_den[:,0],y=(Num_dens - 10**target_num_den[:,1]) / 10**target_num_den[:,1] ,yerr = 0,color="C0")
    axs[int(np.floor(i/2)),i%2].set_xlabel("Magnitude")
    axs[int(np.floor(i/2)),i%2].set_ylabel("Frac Diff")
    axs[int(np.floor(i/2)),i%2].legend(title="Number Density")
    
    fig.savefig(save_path+"_params.png",bbox_inches="tight")
    
    return 0

def plot_HOD_params_errors(HOD_params):
    """
    Plot HOD parameters against magnitude
    """

    labels = ["Mmin","sigma_logm","M0","M1","alpha"]
    fig, axs = plt.subplots(ncols=2,nrows=3,figsize=(7,10))
    for i in range(5):
        axs[int(np.floor(i/2)),i%2].errorbar(x=-target_num_den[:,0],y=HOD_params[:,i],yerr = 0,color="C0")
        axs[int(np.floor(i/2)),i%2].set_xlabel("Magnitude")
        axs[int(np.floor(i/2)),i%2].set_ylabel("Parameter Value")
        axs[int(np.floor(i/2)),i%2].legend(title=labels[i])
    i = 5
    axs[int(np.floor(i/2)),i%2].errorbar(x=-target_num_den[:,0],y=(Num_dens - 10**target_num_den[:,1]) / 10**target_num_den[:,1] ,yerr = 0.0625 * (0.01**0.5) * 1.5 * 3,color="C0")
    axs[int(np.floor(i/2)),i%2].set_xlabel("Magnitude")
    axs[int(np.floor(i/2)),i%2].set_ylabel("Frac Diff")
    axs[int(np.floor(i/2)),i%2].legend(title="Number Density")

    fig.savefig(save_path+"_params_errors.png",bbox_inches="tight")

    return 0


def max_like_params(sampler):
    """Get the parameters which provide the maximum likelihood from the sampling"""

    # Only take after the first 1000 steps to allow burn in
    flat_samples = sampler.backend.get_chain()[:,:,:]
    likelihoods = sampler.backend.get_log_prob()[:,:]

    # Find the maximum likelihood parameter position

    print(np.shape(likelihoods))
    best_param_index1 = np.argmax(likelihoods,axis = 0)
    best_param_index2 = np.argmax(samplers[0].get_log_prob()[best_param_index1,np.arange(num_walkers)])
    best_params = flat_samples[best_param_index1[best_param_index2],best_param_index2,:]
    print(best_param_index1[best_param_index2],best_param_index2)
    print(samplers[0].get_log_prob()[best_param_index1[best_param_index2],best_param_index2])
    return best_params


def max_like_log_prob(sampler):
    """Get the log_probability of the maximum likelihood parameters
    """
    likelihoods = sampler.backend.get_log_prob()[:,:]
    max_like = np.max(likelihoods)
    return max_like
    
    

# Run this directly to perform the fitting:

if __name__ == "__main__":

    # First get the number of halos and the hmf and the fine grained mass bins

    mass_bins_big, cen_halos_big, sat_halos_big = calc_hmf_more_files(path, num_mass_bins_big, mass_bin_edges)

    np.save("mass_bins_big.npy",mass_bins_big)
    np.save("cen_halos_big.npy",cen_halos_big)
    np.save("sat_halos_big.npy",sat_halos_big)
    #mass_bins_big = np.load("mass_bins_big.npy")
    #cen_halos_big = np.load("cen_halos_big.npy")
    #sat_halos_big = np.load("sat_halos_big.npy")

    # Now get the finer grained mass bin centres for accuarate HOD estimation

    mass_bin_centres_big = mass_bins_big[:-1] + np.diff(mass_bins_big)/2
    
    # Now load in the mass - r binned paircounts:

    cencen = np.load(run_path + "_cencen.npy")
    censat = np.load(run_path + "_censat.npy")
    satsat = np.load(run_path + "_satsat.npy")
    satsat_onehalo = np.load(run_path + "_satsat_onehalo.npy")
    
    # Undo the subsampling

    cencen = undo_subsampling(cencen,subsample_array)
    censat = undo_subsampling(censat,subsample_array)
    satsat = undo_subsampling(satsat,subsample_array)
    # Set the random seed

    np.random.seed(random_seed)
    
    print(create_hod_params(initial_params))
    print(priors)
    
    # Do the fitting for each magnitude limit:
    samplers = []
    for i in range(1):
        samples = perform_fitting(i)
        samplers.append(samples)
        
        # Print the parameters at the end of the chain to check they are reasonable
    #    print("Parameters at the end of the fitting chain: (These aren't best fits, just a sanity check)")
    #    print(samples.backend.get_chain()[-1,0,:])

    print("All done!")
    print("Saving outputs")

    # Take the parameters and probabilities and save these as numpy arrays
    # Alternative is pickling the samplers object but can be hard to unpickle
    # if the modules used aren't identical
    
    samplers_to_save_params = np.zeros((num_steps,num_walkers,np.shape(initial_params)[0],1))
    samplers_to_save_probs = np.zeros((num_steps,num_walkers,1))
    for i in range(1):
        samplers_to_save_params[:,:,:,i] = samplers[i].backend.get_chain()
        samplers_to_save_probs[:,:,i] = samplers[i].backend.get_log_prob()
    
    np.save(save_path+"_params.npy",samplers_to_save_params)
    np.save(save_path+"_log_probs.npy",samplers_to_save_probs)


    # Plot the results
    # Find the max likelihood parameters

    best_params = np.zeros((1,len(initial_params)))
    for i in range(1):
        best_params[i,:] = max_like_params(samplers[i])
        print(best_params[i,:])
    
    best_log_prob = max_like_log_prob(samplers[i])
    
    # Plot the HODs and Correlation Function Ratios - 
    # plotting parameters is hard without knowledge of how many there are
    
    HOD_params = create_hod_params(best_params[0,:])
    
    # HOD_params = np.genfromtxt(save_path+"_best_params.txt")
    HODs = np.zeros((len(mass_bin_edges)-1,11))
    CF_ratio = np.zeros((len(r_bin_edges)-1,9))
    CFs = np.zeros((len(r_bin_edges)-1,9))    
    Num_dens = np.zeros(11)
    for i in range(9):
        params = HOD_params[i,:]
        hod_cen_big = Cen_HOD(params,mass_bin_centres_big)
        hod_sat_big = Sat_HOD(params,hod_cen_big,mass_bin_centres_big)


        hod_cen = create_accurate_HOD(hod_cen_big,cen_halos_big,mass_bin_edges,num_mass_bins_big)
        hod_sat = create_accurate_HOD(hod_sat_big,sat_halos_big,mass_bin_edges,num_mass_bins_big)
        cf = create_2pcf(hod_cen,hod_sat,cencen,censat,satsat,satsat_onehalo,num_sat_parts,
                hod_cen_big,hod_sat_big,cen_halos_big,sat_halos_big,
                boxsize)
        HODs[:,i] = hod_cen + hod_sat
        CF_ratio[:,i] = cf / target_2pcf[:,i+1]
        CFs[:,i] = cf
        #print(CF_ratio)

        num_den = calc_number_density(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize)
        Num_dens[i] = num_den

    for i in range(9,11):
        params = HOD_params[i,:]
        hod_cen_big = Cen_HOD(params,mass_bin_centres_big)
        hod_sat_big = Sat_HOD(params,hod_cen_big,mass_bin_centres_big)


        hod_cen = create_accurate_HOD(hod_cen_big,cen_halos_big,mass_bin_edges,num_mass_bins_big)
        hod_sat = create_accurate_HOD(hod_sat_big,sat_halos_big,mass_bin_edges,num_mass_bins_big)
        HODs[:,i] = hod_cen + hod_sat
        num_den = calc_number_density(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize)
        Num_dens[i] = num_den

    r_bin_centres = 10**(np.log10(r_bin_edges[:-1]) + np.diff(np.log10(r_bin_edges))/2 )
    cf_to_save = np.zeros((len(r_bin_centres),10))
    cf_to_save[:,0] = r_bin_centres
    cf_to_save[:,1:] = CFs
    np.savetxt(save_path+"_CF.txt",cf_to_save)
    np.savetxt(save_path+"_num_den.txt",Num_dens)
    plot_CFs(CF_ratio)
    plot_CFs_errors(CF_ratio)
    plot_HODs(HODs)
    plot_HOD_params(HOD_params)
    plot_HOD_params_errors(HOD_params)
    np.savetxt(save_path+"_best_params.txt",HOD_params)
    np.savetxt(save_path+"_best_params_smooth_functions.txt",best_params[0,:])

