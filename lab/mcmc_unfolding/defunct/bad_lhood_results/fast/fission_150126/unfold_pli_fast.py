""" unfolding fast end of proton-lithium spectrum
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")

# path to files json containing paths for all the data
files_json = 'files.json'

# how many parameters in the model
nparam = 9

# names of parameters
param_names = ['mean1','sigma1','peak1',
               'mean2','sigma2','peak2',
               'alpha3','beta3','peak3']

# set initial guesses for each parameter
guesses = [13.5, 0.25, 3e5, 
           9, 0.8, 1e5,
           0.52, 1.1, 5e5]

# load the unfolding data and parameter info
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

# crop the input data if needed
# gs: [112]=1 MeV, [70]=0.1 , [47]=0.01 , [98]=0.5
unfold.group_structure = unfold.group_structure[70:-5]
unfold.reaction_rates = unfold.reaction_rates[2:]
unfold.reaction_rate_uncerts = unfold.reaction_rate_uncerts[2:]
unfold.response_matrix = [i[70:-5] for i in unfold.response_matrix[2:]]
unfold.response_matrix_uncerts = [i[70:-5] for i in 
                                  unfold.response_matrix_uncerts[2:]]

# create neutron flux model, given parameters theta and a neutron energy
# should return neutron flux for a given neutron energy/energies
def model(theta,energy):

    # unpack the tuple of parameters
    (mean1,sigma1,peak1,
     mean2,sigma2,peak2,
     alpha3,beta3,peak3) = theta

    # create the model
    flux = (unfold.gaussian(mean1,sigma1,peak1,energy)+ 
            unfold.gaussian(mean2,sigma2,peak2,energy)+
            unfold.exponential(alpha3,beta3,peak3,energy))
    return flux


# create log-prior, given parameters theta and the neutron flux model
# should return log-prior distribution for the entire group structure
def log_prior(theta,model): 

    # unpack the tuple of parameters
    (mean1,sigma1,peak1,
     mean2,sigma2,peak2,
     alpha3,beta3,peak3) = theta
    
    #mode3= (alpha3-1)*beta3
    #mode3 = beta3*((alpha3-1)/beta3)**(1/alpha3)

    # set hard limits for the parameters
    if (13 < mean1 < 14 and 0.1 < sigma1 < 0.3 and 2e5 < peak1 < 1e9 
        and 7 < mean2 < 11 and 0.6 < sigma2 < 2 and 1e3 < peak2 < 1e9
        and 0.5< alpha3 <0.7 and 1.05<beta3<1.18 and 1e3 < peak3 < 1e6

        #spectrum limits
        and model(theta,15)<1e1 and model(theta,1e-6)<1e3):

        #extra limits for lognormal
        #and 1e0<unfold.lognormal(mode3,sigma3,peak3,10)<1e1):
        
        # extra limits for gamma dist
        #and 0.1<mode3<2 ):

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
nwalkers = 100

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
        params,stds = unfold.postpro_sampler(samples,'corner_fast')
        np.savetxt('results_fast.txt', np.transpose([guesses,params,stds]))