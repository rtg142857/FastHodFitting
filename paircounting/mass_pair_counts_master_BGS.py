# Import modules
import numpy as np
import h5py
import sys
import time
import Corrfunc
from Corrfunc.theory.DD import DD
import fasthod
import config
print("reading in data")

# Load parameters from config file
path = config.path
boxsize = config.boxsize
r_bin_edges = config.r_bin_edges
mass_bin_edges = config.mass_bin_edges
num_sat_parts = config.num_sat_parts
run_label = config.run_label


# In this case we have an hdf5 file so read in using h5py


x, y, z, Mvir, is_central, halo_id = fasthod.read_hdf5(path)

x, y, z, Mvir, x_sat, y_sat, z_sat, Mvir_sat = fasthod.split_cen_sat(x,y,z,Mvir,is_central,halo_id)


import multiprocessing

num_threads =  multiprocessing.cpu_count()
print("CPUS total ",num_threads)
print("total particles ",len(x))

start_time = time.time()
samples_test = fasthod.mass_mask(x,y,z,Mvir,mass_bin_edges)
end_time_1 = time.time()
print('starting pair counting')
npairs_test = fasthod.create_npairs_corrfunc(samples_test,samples_test,r_bin_edges,boxsize,num_threads)
npairs_mass_r_bins_test = fasthod.npairs_conversion(samples_test,samples_test,npairs_test,r_bin_edges)
print('pair counting done')
end_time_2 = time.time()
np.save(run_label+"_cencen.npy",npairs_mass_r_bins_test)

print(run_label)
print(path)
print("time to split into mass bins: ",end_time_1 - start_time)
print("total time to run code including paircounting: ",end_time_2 - start_time)
