import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import config
import config_fitting
# Load parameters from config file
mass_bin_edges = config.mass_bin_edges
mass_bin_centres = config.mass_bin_centres
# Load parameters from fitting config file

Cen_HOD = config_fitting.Cen_HOD
spline_kernel_integral = config_fitting.spline_kernel_integral
cumulative_spline_kernel = config_fitting.cumulative_spline_kernel
Sat_HOD = config_fitting.Sat_HOD



path = config.path
boxsize = config.boxsize
r_bin_edges = config.r_bin_edges
mass_bin_edges = config.mass_bin_edges
num_sat_parts = config.num_sat_parts
run_label = config.run_label

run_path = config_fitting.run_path
num_mass_bins_big = config_fitting.num_mass_bins_big

import fasthod_fitting

create_accurate_HOD = fasthod_fitting.create_accurate_HOD

create_2pcf = fasthod_fitting.create_2pcf

create_weighting_factor = fasthod_fitting.create_weighting_factor

create_randoms_for_corrfunc = fasthod_fitting.create_randoms_for_corrfunc

calc_number_density = fasthod_fitting.calc_number_density

undo_subsampling = fasthod_fitting.undo_subsampling
def var_sample(a,b): # slope = a, complete sample mass = b
    subsampling = 10**(a*(np.log10(mass_bin_centres) -  b))
    subsampling[subsampling > 1] = 1
    return subsampling


subsample_array = var_sample(2,12)

def calc_number_density_array(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize):
    """
    A function which caluculates the number density of halos from the hod
    and a mass function
    """
    num_density = (hod_cen_big*cen_halos_big + hod_sat_big*cen_halos_big) / (boxsize**3)

    return num_density

def convert_num_density_array(num_density,mass_bin_edges,num_mass_bins_big):
    """
    Change big number density array into one using wider mass bins
    """
    num_density_small = np.zeros(len(mass_bin_edges)-1)
    for i in range(len(mass_bin_edges)-1):
        num_density_small = np.sum(num_density[i * num_mass_bins_big / (len(mass_bin_edges)-1):(i+1) * num_mass_bins_big / (len(mass_bin_edges)-1)])
    return num_density_small
if __name__=="__main__":
    """
    Small script to check the HODs and paircounting are being implemented correctly
    Input your HOD parameters on the next line and check the plot looks
    as expected
    Also prints the 2pcf for comparison against an expected one for known paramters
    """
    params = np.array([13.69808703,  0.65566971 ,13.65629438, 14.41252024  ,1.44358639])
    
    mass_bins_big = np.load("mass_bins_big.npy")
    cen_halos_big = np.load("cen_halos_big.npy")
    sat_halos_big = np.load("sat_halos_big.npy")

    # Now get the finer grained mass bin centres for accuarate HOD estimation

    mass_bin_centres_big = mass_bins_big[:-1] + np.diff(mass_bins_big) / 2

    
    
    hod_cen_big = Cen_HOD(params,mass_bin_centres_big)
    hod_sat_big = Sat_HOD(params,hod_cen_big,mass_bin_centres_big)

    cen = Cen_HOD(params,mass_bin_edges)
    sat = Sat_HOD(params,cen,mass_bin_edges)

    
    hod_cen_big = np.zeros(18000)
    hod_cen_big[0:3000] = 1
    hod_sat_big = np.zeros(18000)
    cen = np.zeros(120)
    cen[0:20] = 1
    sat = np.zeros(120)
    
    hod_cen = create_accurate_HOD(hod_cen_big,cen_halos_big,mass_bin_edges,num_mass_bins_big)
    hod_sat = create_accurate_HOD(hod_sat_big,sat_halos_big,mass_bin_edges,num_mass_bins_big)
    
    cencen = np.load(run_path + "_cencen.npy")
    censat = np.load(run_path + "_censat.npy")
    satsat = np.load(run_path + "_satsat.npy")
    satsat_onehalo = np.load(run_path + "_satsat_onehalo.npy")
    
    cencen = undo_subsampling(cencen,subsample_array)
    censat = undo_subsampling(censat,subsample_array)
    satsat = undo_subsampling(satsat,subsample_array)
    
    
    cf = create_2pcf(hod_cen,hod_sat,cencen,censat,satsat,satsat_onehalo,num_sat_parts,
                hod_cen_big,hod_sat_big,cen_halos_big,sat_halos_big,
                boxsize)

    number_density = calc_number_density(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize)
    
    print(subsample_array)
    
    np.savetxt("Check_2pcf_M10_11.txt",cf)
    print(cf)
    print(number_density)
    print(hod_cen)
    print(hod_sat)
    print(np.sum(hod_cen_big*cen_halos_big + hod_sat_big*cen_halos_big))
    #plt.figure(figsize=(8,8))
    #plt.plot(mass_bin_edges,cen,label = "Central")
    #plt.plot(mass_bin_edges,sat,label = "Satellite")
    #plt.plot(mass_bin_edges,cen+sat,label = "Total")
    #plt.ylim(bottom = 1e-4) # Low number cutoff
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.xlabel("Halo Mass /Solar Masses")
    #plt.ylabel("Number of galaxies per halo")
    #plt.legend()
    #plt.show()
    #plt.close()
    


