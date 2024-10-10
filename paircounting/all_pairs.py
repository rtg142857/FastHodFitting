"""Runs all four paircounting things (cencen, censat, satsat, satsat_onehalo), so the data only has to be read in once.
Argument: path_config.yml
"""
# Import modules
import numpy as np
import h5py
import sys
import time
import Corrfunc
from Corrfunc.theory.DD import DD
import fasthod
import config
import yaml
import os
print("reading in data", flush=True)

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

# Read in files using h5py
x, y, z, Mvir, is_central, halo_id = fasthod.read_hdf5_more_files(path, wp_flag, Om0=Om0, Ol0=Ol0, boxSize=boxsize, z_snap=z_snap)

x, y, z, Mvir, x_sat_uncut, y_sat_uncut, z_sat_uncut, Mvir_sat_uncut = fasthod.split_cen_sat(x,y,z,Mvir,is_central,halo_id)

# del x_sat
# del y_sat
# del z_sat
# del Mvir_sat
# del halo_id
# del is_central

x_sat = x_sat_uncut[::num_sat_parts]
y_sat = y_sat_uncut[::num_sat_parts]
z_sat = z_sat_uncut[::num_sat_parts]
Mvir_sat = Mvir_sat_uncut[::num_sat_parts]

x_t, y_t, z_t, Mvir_t = fasthod.read_hdf5_more_files_unresolved(path, wp_flag, Om0=Om0, Ol0=Ol0, boxSize=boxsize, z_snap=z_snap)

x = np.append(x,x_t)
y = np.append(y,y_t)
z = np.append(z,z_t)
Mvir = np.append(Mvir,Mvir_t)

# Starting the paircounting

import multiprocessing

num_threads =  multiprocessing.cpu_count()
print("CPUS total ",num_threads)
print("total particles ",len(x))

# Splitting into mass bins
start_time = time.time()
samples_test = fasthod.mass_mask(x,y,z,Mvir,mass_bin_edges)
samples_test = fasthod.subsample(samples_test,subsample_array)

samples_sat = fasthod.mass_mask(x_sat,y_sat,z_sat,Mvir_sat,mass_bin_edges)
samples_sat = fasthod.subsample(samples_sat,subsample_array)
end_time_1 = time.time()
print("time to split into mass bins: "+str(end_time_1 - start_time), flush=True)

################ Cencen ###############

if not os.path.isfile(run_label+"_cencen.npy"):
    print("============= Cencen pair counting ============", flush=True)
    start_time = time.time()

    print('starting pair counting')
    if not wp_flag:
        npairs_test = fasthod.create_npairs_corrfunc(samples_test,samples_test,r_bin_edges,boxsize,num_threads)
        npairs_mass_r_bins_test = fasthod.npairs_conversion(samples_test,samples_test,npairs_test,r_bin_edges)
    else:
        npairs_test = fasthod.create_npairs_corrfunc_wp(samples_test,samples_test,r_bin_edges,boxsize,num_threads,pi_max,d_pi)
        npairs_mass_r_bins_test = fasthod.npairs_conversion_wp(samples_test,samples_test,npairs_test,r_bin_edges,pi_max, d_pi)
    print('pair counting done')
    end_time_2 = time.time()
    np.save(run_label+"_cencen.npy",npairs_mass_r_bins_test)

    print("Saved with label "+run_label)
    print(path)
    print("total time to run code including paircounting: ",end_time_2 - start_time)
else:
    print("Cencen already done, skipping")

################# Censat ################

if not os.path.isfile(run_label+"_censat.npy"):
    print("============= Censat pair counting ============", flush=True)
    start_time = time.time()

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

    print("Saved with label "+run_label)
    print(path)
    print("total time to run code including paircounting: ",end_time_2 - start_time)
else:
    print("Censat already done, skipping")

################# Satsat ################

if not os.path.isfile(run_label+"_satsat.npy"):
    print("============= Satsat pair counting ============", flush=True)
    start_time = time.time()

    print('starting pair counting')
    if not wp_flag:
        npairs_test = fasthod.create_npairs_corrfunc(samples_sat,samples_sat,r_bin_edges,boxsize,num_threads)
        npairs_mass_r_bins_test = fasthod.npairs_conversion(samples_sat,samples_sat,npairs_test,r_bin_edges)
    else:
        npairs_test = fasthod.create_npairs_corrfunc_wp(samples_sat,samples_sat,r_bin_edges,boxsize,num_threads,pi_max,d_pi)
        npairs_mass_r_bins_test = fasthod.npairs_conversion_wp(samples_sat,samples_sat,npairs_test,r_bin_edges,pi_max, d_pi)

    print('pair counting done')
    end_time_2 = time.time()
    np.save(run_label+"_satsat.npy",npairs_mass_r_bins_test)

    print("Saved with label "+run_label)
    print(path)
    print("total time to run code including paircounting: ",end_time_2 - start_time)
else:
    print("Satsat already done, skipping")

################# Satsat_onehalo ################

if not os.path.isfile(run_label+"_satsat_onehalo.npy"):
    print("============= Satsat_onehalo pair counting ============", flush=True)
    start_time = time.time()
    print('starting pair counting')
    if not wp_flag:
        npairs_test = fasthod.npairs_satsat_onehalo(x_sat_uncut,y_sat_uncut,z_sat_uncut,Mvir_sat_uncut,num_sat_parts,mass_bin_edges,r_bin_edges)
    else:
        npairs_test = fasthod.npairs_satsat_onehalo_wp(x_sat_uncut,y_sat_uncut,z_sat_uncut,Mvir_sat_uncut,num_sat_parts,mass_bin_edges,r_bin_edges,pi_max, d_pi)
    print('pair counting done')
    end_time_2 = time.time()
    np.save(run_label+"_satsat_onehalo.npy",npairs_test)

    print("Saved with label "+run_label)
    print(path)
    print("total time to run code including paircounting: ",end_time_2 - start_time)
else:
    print("Satsat_onehalo already done, skipping")