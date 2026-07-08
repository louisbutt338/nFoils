""" unfolding full proton-lithium spectrum
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")

# path to files json containing paths for all the data
files_json = 'files.json'

# how many parameters in the model
nparam = 11

# names of parameters
param_names = ['mean1','sigma1','peak1',
               'mean2','sigma2','peak2',
               'mode3','peak3',
               'alpha4','beta4','scale4']

# set initial guesses for each parameter
guesses = [13.5, 0.25, 5e5, 
           9, 1, 1e5,
           1, 5e5,
           1.1, 1.2, 5e3]

# load the unfolding data and parameter info
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

# crop the input data if needed
# gs: [112]=1 MeV, [70]=0.1 , [47]=0.01 , [98]=0.5, [5]=1e-6
unfold.group_structure = unfold.group_structure[5:]
unfold.reaction_rates = unfold.reaction_rates
unfold.reaction_rate_uncerts = unfold.reaction_rate_uncerts
unfold.response_matrix = [i[5:] for i in unfold.response_matrix]
unfold.response_matrix_uncerts = [i[5:] for i in 
                                  unfold.response_matrix_uncerts]

# create neutron flux model, given parameters theta and a neutron energy
# should return neutron flux for a given neutron energy/energies
def model(theta,energy):

    # unpack the tuple of parameters
    (mean1,sigma1,peak1,
     mean2,sigma2,peak2,
     mode3,peak3,
     alpha4,beta4,scale4) = theta

    # create the model
    flux = (unfold.gaussian(mean1,sigma1,peak1,energy)+ 
            unfold.gaussian(mean2,sigma2,peak2,energy)+
            unfold.evaporation(mode3,peak3,energy)+
            unfold.epithermal(alpha4,beta4,scale4,1e-7,energy))
    return flux

# create log-prior, given parameters theta and the neutron flux model
# should return log-prior distribution for the entire group structure
def log_prior(theta,model): 

    # unpack the tuple of parameters
    (mean1,sigma1,peak1,
     mean2,sigma2,peak2,
     mode3,peak3,
     alpha4,beta4,scale4) = theta

    # set hard limits for the parameters
    if (13 < mean1 < 14 and 0.1 < sigma1 < 0.3 and 2e5 < peak1 < 1e9 
        and 7 < mean2 < 11 and 0.6 < sigma2 < 2 and 1e3 < peak2 < 1e9
        and 0.2< mode3 <1 and 1e3 < peak3 < 1e6
        and 1<alpha4<1.5 and 1<beta4<2 and 1e1<scale4<1e5

        #spectrum limits
        and model(theta,15)<1e1 and model(theta,1e-7)<1e3):

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
nwalkers = 50

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
        params,stds = unfold.postpro_sampler(samples,'corner')
        np.savetxt('results.txt', np.transpose([guesses,params,stds]))