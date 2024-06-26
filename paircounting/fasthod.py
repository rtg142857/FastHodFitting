# fasthod.py
"""
A module for all the functions used in my code to make mass binned
paircounts from halo catalogs
This should be imported at the start of any scripts
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
import sys
import os
import time
import Corrfunc
from Corrfunc.theory.DD import DD
from Corrfunc.theory.DDrppi import DDrppi


def read_hdf5(path):
    """
    Read in the data from the halo catalog in hdf5
    format at the location specified by path.

    Uses the specific hdf5 format provided to me 
    by Alex Smith when fitting the AbacusSummit
    catalogs for a BGS mock

    Uses h5py
    """

    snap = h5py.File(path,"r")
    position = snap["/position"]
    x = position[:,0]
    y = position[:,1]
    z = position[:,2]
    Mvir = snap["/mass"][:]
    is_central = snap["/is_central"][:]
    halo_id = snap["/halo_id"][:]

    return x, y, z, Mvir, is_central, halo_id

def read_hdf5_wp(path):
    """
    Read in the data from the halo catalog in hdf5
    format at the location specified by path.

    Uses the specific hdf5 format provided to me 
    by Alex Smith when fitting the AbacusSummit
    catalogs for a BGS mock

    Uses h5py
    """

    snap = h5py.File(path,"r")
    position = snap["/position"]
    x = position[:,0]
    y = position[:,1]
    z = position[:,2]
    velocity = snap["/velocity"]
    vz = velocity[:,2]
    Mvir = snap["/mass"][:]
    is_central = snap["/is_central"][:]
    halo_id = snap["/halo_id"][:]
    
    # apply RSD along the z-direction
    Hz = 100.0*np.sqrt(Om0*(1.0+z_snap)**3 + Ol0)
    Hzi = (1+z_snap)/Hz
    z += vz * Hzi

    # apply periodic boundary conditions
    s = z < 0
    z[s] += boxsize
    s = z >= boxsize
    z[s] -= boxsize
    
    return x, y, z, Mvir, is_central, halo_id

def read_hdf5_unresolved_tracers(path):
    """
    Read in position data for tracers representing halos
    below the mass resolution limit which are taken
    from random field particles
    """
    snap = h5py.File(path,"r")
    position = snap["/position"]
    x = position[:,0]
    y = position[:,1]
    z = position[:,2]
    Mvir = snap["/mass"][:]
    
    return x, y, z, Mvir

def read_hdf5_more_files(path, wp_flag=False, Om0=None, Ol0=None, boxSize=None, z_snap=None):
    """
    New halo files are now split into 34 separate files, 
    need a new read routine to get them all in.
    Halo ID is not unique between files 
    """
    
    snap = h5py.File(path+"galaxy_tracers_0.hdf5","r")
    position = snap["/position"]
    x = position[:,0]
    y = position[:,1]
    z = position[:,2]
    Mvir = snap["/mass"][:]
    is_central = snap["/is_central"][:]
    halo_id = snap["/halo_id"][:]

    if wp_flag:
        velocity = snap["/velocity"][:]
        vz = velocity[:,2]

        # apply RSD along the z-direction
        Hz = 100.0*np.sqrt(Om0*(1.0+z_snap)**3 + Ol0)
        Hzi = (1+z_snap)/Hz
        z += vz * Hzi

        # apply periodic boundary conditions
        s = z < 0
        z[s] += boxSize
        s = z >= boxSize
        z[s] -= boxSize
    
    ### FOR FLAMINGO, ONLY USE ONE FILE FOR NOW
    if os.path.exists(path+"galaxy_tracers_1.hdf5", "r"):
        raise Exception("Multiple output files aren't yet supported, fix this")

    # # There won't be double the number of particles
    # # in another file compared to the first one
    # # So use this to create unique halo IDs
    # id_unique = 2 * len(halo_id)
    # for i in range(1,34):
    #     snap = h5py.File(path+"galaxy_tracers_"+str(i)+".hdf5","r")
    #     position_temp = snap["/position"]
    #     x_temp = position_temp[:,0]
    #     y_temp = position_temp[:,1]
    #     z_temp = position_temp[:,2]
    #     Mvir_temp = snap["/mass"][:]
    #     is_central_temp = snap["/is_central"][:]
    #     halo_id_temp = snap["/halo_id"][:] + (i * id_unique)
        
    #     if wp_flag:
    #         velocity_temp = snap["/velocity"][:]
    #         vz_temp = velocity_temp[:,2]

    #         # apply RSD along the z-direction
    #         Hz = 100.0*np.sqrt(Om0*(1.0+z_snap)**3 + Ol0)
    #         Hzi = (1+z_snap)/Hz
    #         z_temp += vz_temp * Hzi

    #         # apply periodic boundary conditions
    #         s = z_temp < 0
    #         z_temp[s] += boxSize
    #         s = z_temp >= boxSize
    #         z_temp[s] -= boxSize
        
    #     x = np.append(x,x_temp)
    #     y = np.append(y,y_temp)
    #     z = np.append(z,z_temp)
    #     Mvir = np.append(Mvir,Mvir_temp)
    #     is_central = np.append(is_central,is_central_temp)
    #     halo_id = np.append(halo_id,halo_id_temp)
    #     print("Reading File Number "+str(i))
    return x, y, z, Mvir, is_central, halo_id


def read_hdf5_more_files_unresolved(path, wp_flag=False, Om0=None, Ol0=None, boxSize=None, 
                                    z_snap=None):
    """
    New halo files are now split into 34 separate files,
    need a new read routine to get them all in.
    Halo ID is not unique between files
    """

    snap = h5py.File(path+"galaxy_tracers_unresolved_0.hdf5","r")
    position = snap["/position"]
    x = position[:,0]
    y = position[:,1]
    z = position[:,2]
    Mvir = snap["/mass"][:]
    
    if wp_flag:
        velocity = snap["/velocity"][:]
        vz = velocity[:,2]

        # apply RSD along the z-direction
        Hz = 100.0*np.sqrt(Om0*(1.0+z_snap)**3 + Ol0)
        Hzi = (1+z_snap)/Hz
        z += vz * Hzi

        # apply periodic boundary conditions
        s = z < 0
        z[s] += boxSize
        s = z >= boxSize
        z[s] -= boxSize

    # There won't be double the number of particles
    # in another file compared to the first one
    # So use this to create unique halo IDs
    for i in range(1,34):
        snap = h5py.File(path+"galaxy_tracers_unresolved_"+str(i)+".hdf5","r")
        position_temp = snap["/position"]
        x_temp = position_temp[:,0]
        y_temp = position_temp[:,1]
        z_temp = position_temp[:,2]
        Mvir_temp = snap["/mass"][:]
    
        if wp_flag:
            velocity_temp = snap["/velocity"][:]
            vz_temp = velocity_temp[:,2]

            ## change to redshift space
            Hz = 100.0*np.sqrt(Om0*(1.0+z_snap)**3 + Ol0)
            Hzi = (1+z_snap)/Hz
            z_temp += vz_temp * Hzi

            # apply periodic boundary conditions
            s = z_temp < 0
            z_temp[s] += boxSize
            s = z_temp >= boxSize
            z_temp[s] -= boxSize
            
        x = np.append(x,x_temp)
        y = np.append(y,y_temp)
        z = np.append(z,z_temp)
        Mvir = np.append(Mvir,Mvir_temp)
        print("Reading Unresolved Halo File Number "+str(i))
    return x, y, z, Mvir



def read_ascii(path):
    """
    Read in the data from a halo catalog in ascii
    format. 

    Uses the specific format provided to me for 
    fitting ELG parameters from DESI simulations

    Includes reading velocities and radii because
    here we need to create satellite particles,
    they are not provided in the catalog
    """
    path = sys.argv[1]
    Mvir, Rvir, Rs, x, y, z, vx, vy, vz, Vrms, PID = np.loadtxt(path,usecols=(2,5,6,8,9,10,11,12,13,4,41),unpack=True)
    

    # Only take the central halos as we shall create our own satellite particles later
    mask1 = np.where(PID == -1)
    
    Mvir = Mvir[mask1]
    Rvir = Rvir[mask1]
    Rs = Rs[mask1]
    x = x[mask1]
    y = y[mask1]
    z = z[mask1]
    vx = vx[mask1]
    vy = vy[mask1]
    vz = vz[mask1]
    Vrms = Vrms[mask1]
    return Mvir, Rvir, Rs, x, y, z, vx, vy, vz, Vrms



def split_cen_sat(x,y,z,Mvir,is_central,halo_id):

    # Split into centrals and satellites

    mask1 = np.array(is_central)
    Mvir_sat = Mvir[~mask1]
    x_sat = x[~mask1]
    y_sat = y[~mask1]
    z_sat = z[~mask1]
    halo_id_sat = halo_id[~mask1]

    halo_id_sat_sorted = np.argsort(halo_id_sat)

    # Sort the satellite arrays according to the parent halo id

    Mvir_sat = Mvir_sat[halo_id_sat_sorted[:]]
    x_sat = x_sat[halo_id_sat_sorted[:]]
    y_sat = y_sat[halo_id_sat_sorted[:]]
    z_sat = z_sat[halo_id_sat_sorted[:]]


    Mvir = Mvir[mask1]
    x = x[mask1]
    y = y[mask1]
    z = z[mask1]

    return x, y, z, Mvir, x_sat, y_sat, z_sat, Mvir_sat


def mass_mask(x,y,z,Mvir,mass_bin_edges):
    """
    Takes arrays of the coordinates and masses of haloes and returns a list
    where the list elements are arrays of halo coordinates filtered into the 
    mass bins provided
    """
    samples_ = []
    for i in range(len(mass_bin_edges)-1):
        mass_mask = (np.array(Mvir>=mass_bin_edges[i]) 
                  & np.array(Mvir<mass_bin_edges[i+1]))
        samples_.append(np.vstack((x[mass_mask],y[mass_mask],z[mass_mask])).T)
    return samples_

def subsample(samples,subsample_array):
    """
    Subsamples the halos contained in the samples array according to the 
    factors contained in the subsample_array
    """
    if len(samples)!= len(subsample_array):
        raise ValueError("number of mass bins is not the same between samples and subsample_array")
    samples_new = []
    for i in range(len(samples)):
        num_part = len(samples[i][:,0])
        subsample_ratio = subsample_array[i]
        randoms = np.random.choice(np.arange(num_part),int(num_part*subsample_ratio),replace = False)
        print(i)
        print(len(samples[i][:,0]))
        samples_new.append(samples[i][randoms,:])
        print(len(samples_new[i][:,0]))
    return samples_new
    
    
    
def create_npairs_corrfunc(samples1,samples2,r_bin_edges,boxsize,num_threads):
    """"
    Takes two lists created from mass_mask(), radial bins, a boxsize,
    and a number of threads to use. 

    From this the number of pairs is counted between every combination
    of mass bin and radial bin

    This paircounting is achieved using corrfunc and can be multithreaded 
    with num_threads

    Returns a list which contains all the paircounts and is then reordered
    into an array in a later function: npairs_conversion
    """
    n_pairs = []
    # Iterate over the length of both mass filtered arrays
    for i in range(len(samples1)):
        for j in range(len(samples2)):
            if len(samples1[i])>1 and len(samples2[j])>1:
                n_pairs.append(DD(0, num_threads, r_bin_edges, samples1[i][:,0], 
                                  samples1[i][:,1],samples1[i][:,2],X2=samples2[j][:,0],
                                  Y2=samples2[j][:,1],Z2 = samples2[j][:,2],
                                  periodic=True, verbose=False,
                                  boxsize=boxsize))
                # We only use Corrfunc if both mass bins are populated, otherwise
                # return 0 for this combination
            else:
                n_pairs.append(0)
            print(i,j)
    return(n_pairs)


def create_npairs_corrfunc_wp(samples1,samples2,r_bin_edges,boxsize,num_threads,pi_max,d_pi=1):
    """"
    Takes two lists created from mass_mask(), radial bins, a boxsize,
    and a number of threads to use. 

    From this the number of pairs is counted between every combination
    of mass bin and radial bin

    This paircounting is achieved using corrfunc and can be multithreaded 
    with num_threads

    Returns a list which contains all the paircounts and is then reordered
    into an array in a later function: npairs_conversion
    """
    n_pairs = []
    # Iterate over the length of both mass filtered arrays
    for i in range(len(samples1)):
        for j in range(len(samples2)):
            if len(samples1[i])>=1 and len(samples2[j])>=1:
                n_pairs.append(DDrppi(autocorr=0, nthreads=num_threads, pimax=pi_max, npibins=(pi_max//d_pi), binfile=r_bin_edges,
                         X1=samples1[i][:,0],Y1=samples1[i][:,1],Z1=samples1[i][:,2],X2=samples2[j][:,0],
                         Y2=samples2[j][:,1],Z2 = samples2[j][:,2],periodic=True,verbose=False, boxsize=boxsize))
                # We only use Corrfunc if both mass bins are populated, otherwise
                # return 0 for this combination
            else:
                n_pairs.append(0)
            print(i,j)
    return(n_pairs)


def npairs_conversion(samples1,samples2,n_pairs,r_bin_edges):
    """
    A simple function to convert the flat list of paircounts from
    create_npairs_corrfunc into a 3D array representing the pairs
    counted in the separate mass and radial bins

    The object returned is an array which contains the paircounts
    for the ith mass bin 1, jth mass bin 2 and, kth r bin
    where we take the [i,j,k] element of the output array
    """
    n_pairs_mass_r_bins = np.zeros((len(samples1),len(samples2),len(r_bin_edges)-1))
    for i in range(len(samples1)):
        for j in range(len(samples2)):
            for k in range(len(r_bin_edges)-1):
                if len(samples1[i])>1 and len(samples2[j])>1:
                    n_pairs_mass_r_bins[i,j,k] = n_pairs[len(samples1)*i + j][k][3]
                else:
                    n_pairs_mass_r_bins[i,j,k] = 0
    return n_pairs_mass_r_bins


def npairs_conversion_wp(samples1,samples2,n_pairs,r_bin_edges,pi_max, d_pi=1):
    """
    A simple function to convert the flat list of paircounts from
    create_npairs_corrfunc into a 3D array representing the pairs
    counted in the separate mass and radial bins

    The object returned is an array which contains the paircounts
    for the ith mass bin 1, jth mass bin 2 and, kth r bin
    where we take the [i,j,k] element of the output array
    """
    n_pairs_mass_r_bins_wp = np.zeros((len(samples1),len(samples2),len(r_bin_edges)-1,(pi_max//d_pi)))
    for i in range(len(samples1)):
        for j in range(len(samples2)):
            for k in range(len(r_bin_edges)-1):
                for l in range((pi_max//d_pi)):
                    if len(samples1[i])>=1 and len(samples2[j])>=1:
                        n_pairs_mass_r_bins_wp[i,j,k,l] = n_pairs[len(samples1)*i + j][k*(pi_max//d_pi) + l][4]
                    else:
                        n_pairs_mass_r_bins_wp[i,j,k,l] = 0
    return n_pairs_mass_r_bins_wp


def npairs_satsat_onehalo(x,y,z,Ms,num_sat_parts,mass_bin_edges,r_bin_edges):
    """
    We cannot use corrfunc for the one halo satellite-satellite term.
    This is because it would be too inefficient to split by
    halo id and use corrfunc on every single halo.

    We can use the fact that the satellite data here is always
    ordered according to halo id 

    We know the number of satellite tracer particles so can use 
    this to count the one halo pairs
    """
    # Create empty arrays to hold data
    # For N sat particles we will have num_sat_parts*(num_sat_parts-1)/2 pairs per halo
    Ms_reduced = np.zeros((len(Ms[::num_sat_parts]),
                           int(num_sat_parts*(num_sat_parts-1)/2)))

    distances = np.zeros((len(Ms[::num_sat_parts]),
                          int(num_sat_parts*(num_sat_parts-1)/2)))
    k = 0 
    # For any number of satellite particles can take every combination of ith and 
    # jth satellite particle and find the distance and put them all in a big array
    for i in range(num_sat_parts):
        for j in range(i):
            print(k)
            Ms_reduced[:,k] = Ms[::num_sat_parts]
            distances[:,k] = ((x[i::num_sat_parts]-x[j::num_sat_parts])**2 
                            + (y[i::num_sat_parts]-y[j::num_sat_parts])**2 
                            + (z[i::num_sat_parts]-z[j::num_sat_parts])**2)**0.5
            # This works as all halos have the same number of sat particles so 
            # [::num_sat_parts] is actually looping the specific paircount 
            # (particles i and j) over every halo
            k += 1
            print(i,j)

    # Reshape the masses and distances of the pairs so they can
    # be binned
    Ms_reduced = np.reshape(Ms_reduced,(1,-1))[0]
    distances = np.reshape(distances,(1,-1))[0]

    final_data = np.histogram2d(x=Ms_reduced,y=distances,bins=[mass_bin_edges,r_bin_edges])
    final_data = final_data[0]
    # Finally transform into the usual format with 2 separate M bins so that it easily fits into the rest of my existing code

    n_pairs_mass_r_bins = np.zeros((len(mass_bin_edges)-1,len(mass_bin_edges)-1,len(r_bin_edges)-1))
    for i in range(len(mass_bin_edges)-1):
        for j in range(len(r_bin_edges)-1):
            n_pairs_mass_r_bins[i,i,j] = final_data[i,j]


    return n_pairs_mass_r_bins


def npairs_satsat_onehalo_wp(x,y,z,Ms,num_sat_parts,mass_bin_edges,r_bin_edges,pi_max, d_pi=1):
    """
    We cannot use corrfunc for the one halo satellite-satellite term.
    This is because it would be too inefficient to split by
    halo id and use corrfunc on every single halo.

    We can use the fact that the satellite data here is always
    ordered according to halo id 

    We know the number of satellite tracer particles so can use 
    this to count the one halo pairs
    """
    # Create empty arrays to hold data
    # For N sat particles we will have num_sat_parts*(num_sat_parts-1)/2 pairs per halo
    Ms_reduced = np.zeros((len(Ms[::num_sat_parts]),
                           int(num_sat_parts*(num_sat_parts-1)/2)))

    # distances in LOS and projected directions are binned
    distances_rp = np.zeros((len(Ms[::num_sat_parts]),
                          int(num_sat_parts*(num_sat_parts-1)/2)))

    distances_pi = np.zeros((len(Ms[::num_sat_parts]),
                          int(num_sat_parts*(num_sat_parts-1)/2)))

    k = 0
    # For any number of satellite particles can take every combination of ith and 
    # jth satellite particle and find the distance and put them all in a big array
    for i in range(num_sat_parts):
        for j in range(i):
            print(k)
            Ms_reduced[:,k] = Ms[::num_sat_parts]
            distances_rp[:,k] = ((x[i::num_sat_parts]-x[j::num_sat_parts])**2
                            + (y[i::num_sat_parts]-y[j::num_sat_parts])**2)**0.5

            distances_pi[:,k] = ((z[i::num_sat_parts]-z[j::num_sat_parts])**2)**0.5

            # This works as all halos have the same number of sat particles so 
            # [::num_sat_parts] is actually looping the specific paircount 
            # (particles i and j) over every halo
            k += 1
            print(i,j)

    # Reshape the masses and distances of the pairs so they can
    # be binned
    Ms_reduced = np.reshape(Ms_reduced,(1,-1))[0]
    distances_rp = np.reshape(distances_rp,(1,-1))[0]
    distances_pi = np.reshape(distances_pi,(1,-1))[0]

    final_data = np.histogramdd(sample = np.array([Ms_reduced,distances_rp,distances_pi]).T,bins=[mass_bin_edges,r_bin_edges,np.arange(0, pi_max+1, d_pi)])
    final_data = final_data[0]
    # Finally transform into the usual format with 2 separate M bins so that it easily fits into the rest of my existing code

    n_pairs_mass_r_bins = np.zeros((len(mass_bin_edges)-1,len(mass_bin_edges)-1,len(r_bin_edges)-1,(pi_max//d_pi)))
    for i in range(len(mass_bin_edges)-1):
        for j in range(len(r_bin_edges)-1):
            for k in range(pi_max//d_pi):
                n_pairs_mass_r_bins[i,i,j,k] = final_data[i,j,k]


    return n_pairs_mass_r_bins


def get_satellites(x,y,z,vx,vy,vz,Vrms,Mvir,Rvir,Rs,boxsize,Ns):
    """
    A function to create Ns satellite particles per halo
    Mostly copy-pasted from Cesar's code

    Places satellite particles randomly according to 
    an NFW profile

    Takes all the arrays of central halo attributes
    as inputs, along with the boxsize and Ns

    Outputs x,y,z,M for the satellite particles
    output arrays are grouped according to the
    parent halo, but in no particular order of
    parent halo.
    """
    c = Rvir/Rs
    fc = np.log(1. + c) - c / (1. + c)
    Rho_s = Mvir / (4*np.pi*Rs**3 * fc)
    Mnfw_r200 = 4. * np.pi * Rho_s * Rs**3 * fc

    x_s = []
    y_s = []
    z_s = []
    vx_s = []
    vy_s = []
    vz_s = []
    Mh_s = []
    Rvir_s = []
    Rs_s = []

    for i in range(len(x)):
        if Ns != 0:
            r   = 10.0**(np.arange(-5.0, 0.05, 0.05)) * Rvir[i]
            theta  = np.arccos(2.0 * np.random.ranf(Ns) - 1.0)
            phi = 2.0 * np.pi * np.random.ranf(Ns)
            M_rand   = Mvir[i] * np.random.ranf(Ns)
            Mnfw = 4. * np.pi * Rho_s[i] * Rs[i]**3 * (np.log(1.+r/Rs[i]) - (r/Rs[i])/(1.+r/Rs[i]))
            rp = np.interp(M_rand, Mnfw, r)
            xp = rp * np.sin(theta) * np.cos(phi)
            yp = rp * np.sin(theta) * np.sin(phi)
            zp = rp * np.cos(theta)

            xs = x[i] + xp/1000. #####Change depending on rockstar catalogue as units can change!!!!
            ys = y[i] + yp/1000. ##### This assumes halo sizes are in kpc/h
            zs = z[i] + zp/1000.

            x_s.append(xs)
            y_s.append(ys)
            z_s.append(zs)

            xx = np.random.normal(0,1,len(xs))
            yy = np.random.normal(0,1,len(ys))
            zz = np.random.normal(0,1,len(zs))

            vxs = vx[i]+ xx*Vrms[i]/np.sqrt(3.)
            vys = vy[i]+ yy*Vrms[i]/np.sqrt(3.)
            vzs = vz[i]+ zz*Vrms[i]/np.sqrt(3.)

            vx_s.append(vxs)
            vy_s.append(vys)
            vz_s.append(vzs)

            Ms = Mvir[i] + (xp - xp)
            Rvirs = Rvir[i] + (xp - xp)
            rs_s = Rs[i] + (xp - xp)

            Mh_s.append(Ms)
            Rvir_s.append(Rvirs)
            Rs_s.append(rs_s)
            #if i%int(len(x)/100)==0:
            #    print(float(i)/len(x))

    x_s = np.hstack(x_s)
    y_s = np.hstack(y_s)
    z_s = np.hstack(z_s)

    vx_s = np.hstack(vx_s)
    vy_s = np.hstack(vy_s)
    vz_s = np.hstack(vz_s)

    Mh_s = np.hstack(Mh_s)
    Rvir_s = np.hstack(Rvir_s)
    Rs_s = np.hstack(Rs_s)

    x_s = x_s.tolist()
    y_s = y_s.tolist()
    z_s = z_s.tolist()

    vx_s = vx_s.tolist()
    vy_s = vy_s.tolist()
    vz_s = vz_s.tolist()

    Mh_s = Mh_s.tolist()
    Rvir_s = Rvir_s.tolist()
    Rs_s = Rs_s.tolist()

    x_s = np.asarray(x_s)
    y_s = np.asarray(y_s)
    z_s = np.asarray(z_s)

    vx_s = np.asarray(vx_s)
    vy_s = np.asarray(vy_s)
    vz_s = np.asarray(vz_s)

    Mh_s = np.asarray(Mh_s)
    Rvir_s = np.asarray(Rvir_s)
    Rs_s = np.asarray(Rs_s)

    for i in range(len(x_s)):
        if x_s[i] < 0:
            x_s[i] = x_s[i] + boxsize
        elif y_s[i] < 0:
            y_s[i] = y_s[i] + boxsize
        elif z_s[i] < 0:
            z_s[i] = z_s[i] + boxsize
        elif x_s[i] > boxsize:
            x_s[i] = x_s[i] - boxsize
        elif y_s[i] > boxsize:
            y_s[i] = y_s[i] - boxsize
        elif z_s[i] > boxsize:
            z_s[i] = z_s[i] - boxsize

    return x_s,y_s,z_s,Mh_s


