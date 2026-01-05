""" example for unfolding spectrum from reaction rates and response functions
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np

# path to files json containing paths for all the data
files_json_path = 'files.json'

# initialise object to load data
unfold = BayesianUnfolding(files_json_path)

# how many parameters in the model
nparam = 3

# names of parameters
param_names = ['mu','sigma','peak']

# get the group structure in MeV so we can create our prior
group_structure = unfold.group_structure

# create neutron flux model, given parameters theta and a neutron energy
def model(self,theta,energy):

    # unpack the tuple of parameters
    mean,sigma,peak = theta

    # create the model
    flux = unfold.gaussian(mean,sigma,peak,energy)
    return flux

# create neutron spectrum prior, given parameters theta
# using the model created above
def prior(self,theta): 

    # unpack the tuple of parameters
    mean,sigma,peak = theta

    # set hard limits for the parameters
    if (12 < mean < 15 and 0.1 < sigma < 2 and 1e3 < peak < 1e9):

        # return model for the entire group structure if inside limits, 
        # or return -inf if outside limits 
        return model(self,theta,group_structure)
    return -np.inf

# bind the new model and prior to the object methods
unfold.model = model.__get__(unfold,BayesianUnfolding)
unfold.prior = prior.__get__(unfold,BayesianUnfolding)

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

# set initial guesses for each parameter
guesses = [13.5, 1.5, 1e7]

# run using the above inputs and plot
params,stds = unfold.run_ensemble_mcmc(nparam,param_names,rm_samples,
                                       nwalkers,nburn,nsteps,guesses)
unfold.plot_spectrum(guesses,params,stds)

