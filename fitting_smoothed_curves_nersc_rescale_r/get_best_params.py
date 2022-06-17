import numpy as np
import sys

samples_in = sys.argv[1]
likelihoods_in = sys.argv[2]
outname = sys.argv[3]

flat_samples = np.load(samples_in)
likelihoods = np.load(likelihoods_in)[:,:,0]

num_walkers = 180
def max_like_params(flat_samples,likelihoods):
    """Get the parameters which provide the maximum likelihood from the sampling"""

    # Find the maximum likelihood parameter position

    print(np.shape(likelihoods))
    best_param_index1 = np.argmax(likelihoods,axis = 0)
    print(best_param_index1)
    best_param_index2 = np.argmax(likelihoods[best_param_index1,np.arange(num_walkers)])
    print(best_param_index2)
    best_params = flat_samples[best_param_index1[best_param_index2],best_param_index2,:]
    print(best_param_index1[best_param_index2],best_param_index2)
    print(likelihoods[best_param_index1[best_param_index2],best_param_index2])
    return best_params


best_params = max_like_params(flat_samples,likelihoods)

np.savetxt(outname,best_params)
print("Done!")
