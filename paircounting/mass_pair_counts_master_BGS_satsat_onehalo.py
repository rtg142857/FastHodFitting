# Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import halotools
import halotools.mock_observables
import h5py
import sys
import time
import Corrfunc
from Corrfunc.theory.DD import DD
import fasthod
import config
print("reading in data")

# Load data from config file
path = config.path
boxsize = config.boxsize
r_bin_edges = config.r_bin_edges
mass_bin_edges = config.mass_bin_edges
num_sat_parts = int(config.num_sat_parts)
run_label = config.run_label

# In this case we have an hdf5 file so read in using h5py


x, y, z, Mvir, is_central, halo_id = fasthod.read_hdf5(path)

x, y, z, Mvir, x_sat, y_sat, z_sat, Mvir_sat = fasthod.split_cen_sat(x,y,z,Mvir,is_central,halo_id)

# Now we only want to take 1 satellite particle per halo
# Currently 3 satellite particles per halo


import multiprocessing

num_threads = 2  # multiprocessing.cpu_count()
print("CPUS total ",num_threads)
print("total particles ",len(x))

start_time = time.time()
print('starting pair counting')
npairs_test = fasthod.npairs_satsat_onehalo(x_sat,y_sat,z_sat,Mvir_sat,num_sat_parts,mass_bin_edges,r_bin_edges)
print('pair counting done')
end_time_2 = time.time()
np.save(run_label+"_satsat_onehalo.npy",npairs_test)

print(run_label)
print(path)
print("total time to run code including paircounting: ",end_time_2 - start_time)
