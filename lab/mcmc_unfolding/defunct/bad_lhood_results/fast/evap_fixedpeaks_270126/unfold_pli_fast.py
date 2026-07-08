""" unfolding fast end of proton-lithium spectrum
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np
import multiprocessing as mp

# path to files json containing paths for all the data
files_json = 'files.json'

# how many parameters in the model
nparam = 8

# names of parameters
param_names = ['sigma1','peak1',
               'sigma2','peak2',
               'sigma3','peak3',
               'mode4','peak4']

# set initial guesses for each parameter
guesses = [0.3, 5e5, 
           0.8, 1e5,
           0.8, 1e5,
           1, 5e5]

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
    (sigma1,peak1,
     sigma2,peak2,
     sigma3,peak3,
     mode4,peak4) = theta

    # create the model
    flux = (unfold.gaussian(15.5-2.0945,sigma1,peak1,energy)+ 
            unfold.gaussian(15.5-6.45,sigma2,peak2,energy)+
            unfold.gaussian(15.5-8.85,sigma3,peak3,energy)+
            unfold.weisskopf(mode4,peak4,energy))
    return flux

# create log-prior, given parameters theta and the neutron flux model
# should return log-prior distribution for the entire group structure
def log_prior(theta,model): 

    # unpack the tuple of parameters
    (sigma1,peak1,
     sigma2,peak2,
     sigma3,peak3,
     mode4,peak4) = theta
    
    #mode3= (alpha3-1)*beta3
    #mode3 = beta3*((alpha3-1)/beta3)**(1/alpha3)

    # set hard limits for the parameters
    if ( 0.1 < sigma1 < 1 and 1e3 < peak1 < 1e9 
         and 0.1 < sigma2 < 1 and 1e3 < peak2 < 1e9
         and 0.1 < sigma3<1 and 1e1<peak3<1e9
         and 0.1< mode4 < 2 and 1e3 < peak4 < 1e6

        #spectrum limits
        and model(theta,15)<1e1): # and model(theta,1e-2)<1e4):
        
        # generate and sum model values to get prior
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
nwalkers = 50

# "burn-in" period to let chains stabilize
nburn = 500

# total number of MCMC steps to take (including nburn)
# no. of trace results =  nwalkers * (nsteps-nburn)
nsteps = 1000

# run sampler, postprocess, save results (single)
#samples = unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
#                             rm_samples,nwalkers,nburn,nsteps)
#params,stds = unfold.postpro_sampler(samples,'corner_test')
#print('BIC =',unfold.get_bic(model,params,stds))
#np.savetxt('results_test.txt', np.transpose([guesses,params,stds]))

# run sampler, postprocess, save results (parallel)
if __name__ == '__main__':
   with mp.Pool() as pool:
        samples=unfold.run_sampler(log_posterior,model,log_prior,log_likelihood,
                                   rm_samples,nwalkers,nburn,nsteps,pool)
        params,stds = unfold.postpro_sampler(samples,'corner_fast2')
        print('BIC =',unfold.get_bic(model,params,stds))
        np.savetxt('results_fast2.txt', np.transpose([guesses,params,stds]))
