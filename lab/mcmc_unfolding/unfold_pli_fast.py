""" unfolding fast end of proton-lithium spectrum
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np
import multiprocessing as mp

# path to files json containing paths for all the data
files_json = 'files.json'

# how many parameters in the model
nparam = 4

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

# crop the input data when using 175 group
# gs: [112]=1 MeV, [70]=0.1 , [47]=0.01 , [98]=0.5
unfold.group_structure = unfold.group_structure[:]
unfold.reaction_rates = unfold.reaction_rates[:]
unfold.reaction_rate_uncerts = unfold.reaction_rate_uncerts[:]
unfold.response_matrix = [i[:] for i in unfold.response_matrix[:]]
unfold.response_matrix_uncerts = [i[:] for i in 
                                  unfold.response_matrix_uncerts[:]]

# adjust input data using linear groups
unfold.response_matrix_uncerts = [i/100 for i in unfold.response_matrix_uncerts]

# create neutron flux model, given parameters theta and a neutron energy
# should return neutron flux for a given neutron energy/energies
def model(theta,energy):

    # unpack the tuple of parameters
    (#mean1,
     #sigma1,peak1,
     #sigma2,peak2,
     #sigma3,peak3,
     mode4,peak4
     ) = theta

    # create the model
    flux = (#unfold.gaussian(mean1,sigma1,peak1,energy)
            #unfold.gaussian(15.5-6.45,sigma2,peak2,energy)+
            #unfold.gaussian(15.5-8.85,sigma3,peak3,energy)+

            # REMEMBER THIS ONE
            unfold.weisskopf(mode4,peak4,energy)
            )
    return flux

# create log-prior, given parameters theta and the neutron flux model
# should return log-prior distribution for the entire group structure
def log_prior(theta,model): 

    # check physical spectrum limits
    if model(theta,16)<1e1:

        # then find the joint uniform probability distribution
        # this should include all the parameter limits
        #  P(param) = 1/(upper_limit-lower_limit)
        # one column for each parameter
        prior = 1/( 
            #(14-13) # mean1
            #*(0.8-0.1) #sigma1
            #* (1e6-5e4) #peak1
            #* (0.8-0.1) #sigma2
            #* (5e5-5e4) #peak2
             (2-0.1) #mode4
            * (5e6-1e5) #peak4
            )
        return np.log(prior)
    else:
        return -np.inf

# grab the log likelihood and log posterior
# these are automatically constructed from the prior/model at runtime
log_likelihood = unfold.log_likelihood
log_posterior = unfold.log_posterior

# number of samples from the response function distributions
# computational expense skyrockets here, so start with 1
rm_samples = 50

# number of MCMC walkers/chains
nwalkers = 50

# "burn-in" period to let chains stabilize
nburn = 50

# total number of MCMC steps to take (including nburn)
# no. of trace results =  rm_samples * nwalkers * (nsteps-nburn)
nsteps = 100

# run sampler, postprocess, save results (single)
#samples = unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
#                             rm_samples,nwalkers,nburn,nsteps)
#np.savetxt('samples.txt',samples)

# run sampler, postprocess, save results (parallel)
if __name__ == '__main__':
   with mp.Pool() as pool:
        samples=unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
                                   rm_samples,nwalkers,nburn,nsteps,pool)
        np.savetxt('samples.txt',samples)
