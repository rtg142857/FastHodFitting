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

target_2pcf = config_fitting.target_2pcf

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
    params = np.genfromtxt("number_density/params_fixed_number_density.txt")
    # params = np.array([13.69808703,  0.65566971 ,13.65629438, 14.41252024  ,1.44358639])
    

    mass_bins_big = np.load("mass_bins_big.npy")
    cen_halos_big = np.load("cen_halos_big.npy")
    sat_halos_big = np.load("sat_halos_big.npy")

    # Now get the finer grained mass bin centres for accuarate HOD estimation

    mass_bin_centres_big = mass_bins_big[:-1] + np.diff(mass_bins_big) / 2

    cencen = np.load(run_path + "_cencen.npy")
    censat = np.load(run_path + "_censat.npy")
    satsat = np.load(run_path + "_satsat.npy")
    satsat_onehalo = np.load(run_path + "_satsat_onehalo.npy")
    number = np.zeros(9)
    plt.figure(figsize=(8,8))
    for i in range(9):
        hod_cen_big = Cen_HOD(params[i,1:],mass_bin_centres_big,None,None)
        hod_sat_big = Sat_HOD(params[i,1:],hod_cen_big,mass_bin_centres_big,None,None,None)

        hod_cen = create_accurate_HOD(hod_cen_big,cen_halos_big,mass_bin_edges,num_mass_bins_big)
        hod_sat = create_accurate_HOD(hod_sat_big,sat_halos_big,mass_bin_edges,num_mass_bins_big)
    
    
        cf = create_2pcf(hod_cen,hod_sat,cencen,censat,satsat,satsat_onehalo,num_sat_parts,
                hod_cen_big,hod_sat_big,cen_halos_big,sat_halos_big,
                boxsize)

        number_density = calc_number_density(hod_cen_big,hod_sat_big,
                        cen_halos_big,boxsize)
        number[i] = number_density
        r_bin_centres = 10**(np.log10(r_bin_edges[:-1]) + np.diff(np.log10(r_bin_edges)) / 2)
        plt.plot(r_bin_centres,cf/target_2pcf[:,i+1],label = params[i,0])
        
    plt.xscale("log")
    # plt.yscale("log")
    plt.ylim(0.5,2)
    plt.grid()
    plt.xlabel("r [Mpch$^{-1}$]")
    plt.ylabel(r"$\xi(r)$")
    plt.legend()
    plt.show()
    plt.close()



    params = np.genfromtxt("number_density/params_fixed_number_density.txt")

    number_target = 10**config_fitting.target_num_den[:,1]
    num_den_diff = (number - number_target) / number_target
    labels = ["Mmin","sigma_logm","M0","M1","alpha"]
    fig, axs = plt.subplots(ncols=2,nrows=3,figsize=(15,20))
    log_params = True
    for i in range(5):
        axs[int(np.floor(i/2)),i%2].plot(params[:,0],params[:,i+1],color="C0")
        if not log_params:
            axs[int(np.floor(i/2)),i%2].set_yscale("log")
        axs[int(np.floor(i/2)),i%2].set_xlabel("Magnitude")
        axs[int(np.floor(i/2)),i%2].set_ylabel("Parameter Value")
        axs[int(np.floor(i/2)),i%2].legend(title=labels[i])
    axs[2,1].plot(params[:,0],num_den_diff)
    axs[2,1].set_xlabel("Magnitude")
    #axs[2,1].set_ylabel("")
    axs[2,1].legend(title="Number density fractional difference")

    plt.show()

    
    

