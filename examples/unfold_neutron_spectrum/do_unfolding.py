""" example for unfolding spectrum from reaction rates and response functions 
in sequential or in parallel
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore") # ignore RuntimeWarnings

# path to files json containing paths for all the data
files_json = 'files.json'

# how many parameters in the model
nparam = 3

# names of parameters
param_names = ['mu','sigma','peak']

# set initial guesses for each parameter
guesses = [13.5, 1.5, 1e7]

# load the unfolding data and parameter info
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

# create neutron flux model, given parameters theta and a neutron energy
# should return neutron flux for a given neutron energy/energies
def model(theta,energy):

    # unpack the tuple of parameters
    mean,sigma,peak = theta

    # create the model
    flux = unfold.gaussian(mean,sigma,peak,energy)
    return flux

# create log-prior, given parameters theta and the neutron flux model
# should return log-prior distribution for the entire group structure
def log_prior(theta,model): 

    # unpack the tuple of parameters
    mean,sigma,peak = theta

    # set hard limits for the parameters
    if (12 < mean < 15 and 0.1 < sigma < 2 and 1e3 < peak < 1e9):

        # return log(model) for the group structure (MeV) if in limits
        logprior = np.log(model(theta,unfold.group_structure))
        return logprior
    
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
nwalkers = 20

# "burn-in" period to let chains stabilize
nburn = 500

# total number of MCMC steps to take (including nburn)
# no. of trace results =  nwalkers * (nsteps-nburn)
nsteps = 5000  

# run sampler, postprocess, save results, plot spectrum (single)
samples = unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
                             rm_samples,nwalkers,nburn,nsteps)
params,stds = unfold.postpro_sampler(samples)
np.savetxt('results.txt', np.transpose([guesses,params,stds]))
unfold.plot_simple_spectrum(model,params,stds)

# run sampler, postprocess, save results, plot spectrum (parallel)
#if __name__ == '__main__':
#   with mp.Pool() as pool:
#        samples=unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
#                                   rm_samples,nwalkers,nburn,nsteps,pool)
#        params,stds = unfold.postpro_sampler(samples)
#        np.savetxt('results.txt',np.transpose([guesses,params,stds]))
#params = np.loadtxt('results.txt', usecols=(1,))
#stds = np.loadtxt('results.txt', usecols=(2,))
#unfold.plot_simple_spectrum(model,params,stds)

