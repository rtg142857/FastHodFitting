# Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
import yaml
import sys
import time
import Corrfunc
from Corrfunc.theory.DD import DD
import fasthod
import config
print("reading in data")
import time
time.sleep(1)

# Load parameters from FastHodFitting config file
path = config.path
r_bin_edges = config.r_bin_edges
mass_bin_edges = config.mass_bin_edges
num_sat_parts = config.num_sat_parts
run_label = config.run_label
subsample_array = config.subsample_array

wp_flag = config.wp_flag
pi_max = config.pi_max
d_pi = config.d_pi

# Load parameters from HOD_Mock_Pipeline config file

path_config_filename = sys.argv[1]
with open(path_config_filename, "r") as file:
    path_config = yaml.safe_load(file)

z_snap = path_config["Params"]["redshift"]
boxsize = path_config["Params"]["L"]

flamingo_param_file_path = path_config["Paths"]["params_path"]
with open(flamingo_param_file_path, "r") as file:
    flamingo_params = yaml.safe_load(file)

Om0 = flamingo_params["Cosmology"]["Omega_cdm"]
Ol0 = flamingo_params["Cosmology"]["Omega_lambda"]

# In this case we have an hdf5 file so read in using h5py


x, y, z, Mvir, is_central, halo_id = fasthod.read_hdf5_more_files(path, wp_flag, Om0=Om0, Ol0=Ol0, boxSize=boxsize, z_snap=z_snap)

x, y, z, Mvir, x_sat, y_sat, z_sat, Mvir_sat = fasthod.split_cen_sat(x,y,z,Mvir,is_central,halo_id)

x_sat = x_sat[::num_sat_parts]
y_sat = y_sat[::num_sat_parts]
z_sat = z_sat[::num_sat_parts]
Mvir_sat = Mvir_sat[::num_sat_parts]
samples_sat = fasthod.mass_mask(x_sat,y_sat,z_sat,Mvir_sat,mass_bin_edges)
samples_sat = fasthod.subsample(samples_sat,subsample_array)
del x_sat
del y_sat
del z_sat
del Mvir_sat
del halo_id
del is_central

x_t, y_t, z_t, Mvir_t = fasthod.read_hdf5_more_files_unresolved(path, wp_flag, Om0=Om0, Ol0=Ol0, boxSize=boxsize, z_snap=z_snap)

x = np.append(x,x_t)
y = np.append(y,y_t)
z = np.append(z,z_t)
Mvir = np.append(Mvir,Mvir_t)


# Now we only want to take 1 satellite particle per halo
# Currently 3 satellite particles per halo


import multiprocessing

num_threads = multiprocessing.cpu_count()
print("CPUS total ",num_threads)
print("total particles ",len(x))

start_time = time.time()

samples_test = fasthod.mass_mask(x,y,z,Mvir,mass_bin_edges)


samples_test = fasthod.subsample(samples_test,subsample_array)


end_time_1 = time.time()
print('starting pair counting')
if not wp_flag:
    npairs_test = fasthod.create_npairs_corrfunc(samples_test,samples_sat,r_bin_edges,boxsize,num_threads)
    npairs_mass_r_bins_test = fasthod.npairs_conversion(samples_test,samples_sat,npairs_test,r_bin_edges)
else:
    npairs_test = fasthod.create_npairs_corrfunc_wp(samples_test,samples_sat,r_bin_edges,boxsize,num_threads,pi_max,d_pi)
    npairs_mass_r_bins_test = fasthod.npairs_conversion_wp(samples_test,samples_sat,npairs_test,r_bin_edges,pi_max,d_pi)
print('pair counting done')
end_time_2 = time.time()
np.save(run_label+"_censat.npy",npairs_mass_r_bins_test)

print(run_label)
print(path)
print("time to split into mass bins: ",end_time_1 - start_time)
print("total time to run code including paircounting: ",end_time_2 - start_time)
