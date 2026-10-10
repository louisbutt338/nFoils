""" postprocess the samples from the emcee results
"""
import numpy as np
from nfoils.unfold import BayesianUnfolding

# path to files json containing paths for all the data
files_json = 'files.json'

# how many parameters in the model
nparam = 2

# names of parameters
param_names = [#r'$\mu_{fast,1}$', 
               #r'$\sigma_{fast,1}$',r'$C_{fast,1}$',
               #r'$\sigma_{fast,2}$',r'$C_{fast,2}$',
               #r'$\sigma_{fast,3}$',r'$C_{fast,3}$',
               r'$T_{evap}$', r'$C_{evap}$'
               ]


# set initial guesses for each parameter
guesses = [#15.5-2.0945,
           #0.35, 5e5, 
           #0.5, 1e5,
           #0.5, 1e4,
           1, 1e6
           ]

# load the unfolding data and parameter info
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

#load samples and do corner plot
samples = np.loadtxt('samples.txt')
params,cov_matrix = unfold.postpro_sampler(samples,'corner')

# save params and covs for further postprocessing
np.savetxt('params.txt', np.transpose(params))
np.savetxt('covs.txt', np.transpose(cov_matrix))