import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import config
import config_fitting
# Load parameters from config file
mass_bin_edges = config.mass_bin_edges

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

    
    hod_cen = create_accurate_HOD(hod_cen_big,cen_halos_big,mass_bin_edges,num_mass_bins_big)
    hod_sat = create_accurate_HOD(hod_sat_big,sat_halos_big,mass_bin_edges,num_mass_bins_big)
    
    cencen = np.load(run_path + "_cencen.npy")
    censat = np.load(run_path + "_censat.npy")
    satsat = np.load(run_path + "_satsat.npy")
    satsat_onehalo = np.load(run_path + "_satsat_onehalo.npy")
    
    cf = create_2pcf(hod_cen,hod_sat,cencen,censat,satsat,satsat_onehalo,num_sat_parts,
                hod_cen_big,hod_sat_big,cen_halos_big,sat_halos_big,
                boxsize)

    number_density = calc_number_density(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize)
    
    print(cen_halos_big[:100])
    print(cf)
    print(number_density)
    print(hod_cen)
    print(hod_sat)
    print(cen_halos_big[10000])
    print(np.sum(hod_cen_big*cen_halos_big + hod_sat_big*cen_halos_big))
    plt.figure(figsize=(8,8))
    plt.plot(mass_bin_edges,cen,label = "Central")
    plt.plot(mass_bin_edges,sat,label = "Satellite")
    plt.plot(mass_bin_edges,cen+sat,label = "Total")
    plt.ylim(bottom = 1e-4) # Low number cutoff
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Halo Mass /Solar Masses")
    plt.ylabel("Number of galaxies per halo")
    plt.legend()
    plt.show()
    plt.close()
    


