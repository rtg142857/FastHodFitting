Create mass-distance binned paircounts of dark matter haloes

In config.py mass bins are specified by mass_bin_edges and pair-distance bins are specified by r_bin_edges

r_bin_edges should be the same as those used in the target correlation function

path should contain the location of the halo catalogue from which the paircounts are to be measured

subsampling of the halos can be applied by changing subsample_array, this should only be done for low mass haloes.

Each set of pairs (central-central, central-satellite, satellite-satellite, satellite-satellite one halo) is calculated
using its own python script. These can be run independently to generate the tabulated paircounts

The paircounts are saved as numpy objects with the prefix of the filename defined by the run_label string in config.py


