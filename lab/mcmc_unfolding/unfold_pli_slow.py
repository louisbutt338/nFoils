""" unfolding slow end of proton lithium spectrum
"""
from bfoils.unfold import BayesianUnfolding
import numpy as np
import multiprocessing as mp

# path to files json containing paths for all the data
files_json = 'files.json'

# how many parameters in the model
nparam = 2

# names of parameters
param_names = ['alpha','scale']

# set initial guesses for each parameter
guesses = [1.1,1e3]

# load the unfolding data and parameter info
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

# crop the input data if needed
# gs: [112]=1 MeV, [70]=0.1 , [47]=0.01 , [98]=0.5, [5]=1e-6
unfold.group_structure = unfold.group_structure[:47]
unfold.reaction_rates = unfold.reaction_rates[:2]
unfold.reaction_rate_uncerts = unfold.reaction_rate_uncerts[:2]
unfold.response_matrix = [i[:47] for i in unfold.response_matrix[:2]]
unfold.response_matrix_uncerts = [i[:47] for i in 
                                  unfold.response_matrix_uncerts[:2]]

# create neutron flux model, given parameters theta and a neutron energy
# should return neutron flux for a given neutron energy/energies
def model(theta,energy):

    # unpack the tuple of parameters
    (alpha,scale) = theta

    # create the model
    flux = unfold.powerlaw(alpha,scale,energy)
    return flux

# create log-prior, given parameters theta and the neutron flux model
# should return log-prior distribution for the entire group structure
def log_prior(theta,model): 

    # unpack the tuple of parameters
    (alpha,scale) = theta
    
    # set hard limits for the parameters
    if (1<alpha<2 and 1e2<scale<1e4):

        # return log(prior) for the group structure (MeV) if in limits
        prior = np.sum([model(theta,i) for i in unfold.group_structure])
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
nsteps = 1000

# run sampler, postprocess, save results (single)
#samples = unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
#                             rm_samples,nwalkers,nburn,nsteps)
#params,stds = unfold.postpro_sampler(samples)
#np.savetxt('results.txt', np.transpose([guesses,params,stds]))

# run sampler, postprocess, save results (parallel)
if __name__ == '__main__':
   with mp.Pool() as pool:
        samples=unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
                                   rm_samples,nwalkers,nburn,nsteps,pool)
        params,stds = unfold.postpro_sampler(samples,'corner_slow')
        print('BIC =',unfold.get_bic(model,params,stds))
        np.savetxt('results_slow.txt', np.transpose([guesses,params,stds]))