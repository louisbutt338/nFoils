""" example for unfolding a 14 MeV peak spectrum 
from reaction rates and response functions with uncertainties. 
options for running in sequential or in parallel
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np
import multiprocessing as mp

# path to files json containing paths for all the unfolding data
# should include response matrix with uncertainties, group structure
# and reaction rates with uncertainties
files_json = 'files.json'

# how many parameters in the model
nparam = 2

# names of parameters
param_names = ['sigma','peak']

# set initial guesses for each parameter
guesses = [0.5, 1e5]

# load the unfolding data and parameter info
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

# create neutron flux model, given parameters theta 
# and a given neutron energy (MeV)
def model(theta,energy):

    # unpack the tuple of parameters
    sigma,peak = theta

    # set fixed mean energy of the gaussian
    mean_energy = 14

    # create the model and return the flux
    flux = unfold.gaussian(mean_energy,sigma,peak,energy)
    return flux

# create log-prior, given parameters theta and the neutron flux model
# should return log-prior distribution for the entire group structure
def log_prior(theta,model): 

    # unpack the tuple of parameters
    sigma,peak = theta

    # set hard limits for the parameters
    if (0.1 < sigma < 2 and 1e3 < peak < 1e9):

        # generate and sum model values for the group structure
        # to get the prior
        prior = np.sum([model(theta,i) for i in unfold.group_structure])

        # return log(prior) for the group structure (MeV) if in limits
        return np.log(prior)
    
    # return -inf if outside limits 
    return -np.inf

# grab the log likelihood and log posterior
# these are automatically constructed from the prior/model at runtime
log_likelihood = unfold.log_likelihood
log_posterior = unfold.log_posterior

# number of samples from the response function distributions
# computational expense skyrockets here, so start with 1
rm_samples = 1

# number of MCMC walkers/chains
nwalkers = 30

# "burn-in" period to let chains stabilize
nburn = 500

# total number of MCMC steps to take (including nburn)
# no. of trace results =  nwalkers * (nsteps-nburn)
nsteps = 5000

# run sampler, postprocess, save results, plot spectrum (single)
#samples = unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
#                             rm_samples,nwalkers,nburn,nsteps)
#params,stds = unfold.postpro_sampler(samples)
#np.savetxt('results.txt', np.transpose([guesses,params,stds]))
#unfold.plot_simple_spectrum(model,params,stds)

# run sampler, postprocess, save results, plot spectrum (parallel)
if __name__ == '__main__':
   with mp.Pool() as pool:
        samples=unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
                                   rm_samples,nwalkers,nburn,nsteps,pool)
        params,stds = unfold.postpro_sampler(samples)
        np.savetxt('results.txt',np.transpose([guesses,params,stds]))
        unfold.plot_simple_spectrum(model,params,stds)