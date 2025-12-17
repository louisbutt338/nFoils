""" 
module for doing Bayesian unfolding with emcee

D Foreman-Mackey et al, emcee: The MCMC Hammer, 
Instrumentation and Methods for Astrophysics 2012
"""

import emcee
import corner
import warnings
import numpy as np
import csv
import json
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=20)

class BayesianUnfolding:
    """ class to do mcmc unfolding 
    """

    def __init__(self,files_json_path): 
        """ Initialise Bayesian unfolding class, and set all unfolding data as 
        numpy arrays of floats

        Attributes
        ----------
        files_json_path : str
            path to the files data json
        """

        # set files json attribute
        self.files_json_path = files_json_path

        # load the datafiles paths
        with open(self.files_json_path
                  ) as files_json:
            files_paths = json.load(files_json)

        # set energy grid attribute
        energy_grid_path = files_paths["group_structure"]
        self.energy_grid = np.fromfile(energy_grid_path, sep=" ")

        # set reaction rate attributes
        with open(files_paths["reaction_rates"]) as file:
            reaction_rate_csv = csv.reader(file,delimiter=',')
            self.reaction_rates = np.array(list(reaction_rate_csv),dtype=float)
        with open(files_paths["reaction_rate_uncerts"]) as file:
            reaction_uncerts_csv = csv.reader(file,delimiter=',')
            self.reaction_rate_uncerts = np.array(list(reaction_uncerts_csv),dtype=float)
        with open(files_paths["reaction_rate_labels"]) as file:
            reaction_labels_csv = csv.reader(file,delimiter=',')
            self.reaction_rate_labels = np.array(list(reaction_labels_csv))

        # set response matrix attributes
        with open(files_paths["response_matrix"]) as file:
            response_matrix_csv = csv.reader(file,delimiter=',')
            self.response_matrix = np.array(list(response_matrix_csv),dtype=float)
        with open(files_paths["response_matrix_uncerts"]) as file:
            response_uncerts_csv = csv.reader(file,delimiter=',')
            self.response_matrix_uncerts = np.array(list(response_uncerts_csv),dtype=float)

        # adjust for testing just the 14 mev peak, with only the last reaction
        # use [162:] gs values and [6:] rr/rf for just 14 MeV peak
        # use [140:] gs values and [3:] rr/rf for 14 and 8 MeV peak
        self.energy_grid = self.energy_grid[140:-1] * 1e-6
        self.reaction_rates = self.reaction_rates[3:]
        self.reaction_rate_uncerts = self.reaction_rate_uncerts[3:]
        self.response_matrix = [i[140:] for i in self.response_matrix[3:]]
        self.response_matrix_uncerts = [i[140:] for i in self.response_matrix_uncerts[3:]]


    # example probability distributions to construct prior and model
    def gaussian(self,diff,peak,sigma):
        return (( peak / np.sqrt(2 * np.pi * sigma ** 2) )
                * np.exp( - 0.5 * diff ** 2 / sigma ** 2))

    # empty functions for model and prior
    def spectrum_model(self,theta,energy):
        """ model for the spectrum, for calculation with the group structure
        should be theta (the parameters) and energy input here
        """
        print('please create spectrum model first')
        exit()

    def spectrum_prior(self,theta):
        """ prior for the spectrum
        """
        print('please create spectrum prior first')
        exit()

    # likelihood for reaction rates
    def _log_likelihood(self,theta, rr, sigma_rr, response):
        """ gaussian distribution for the reaction rate measurements

        input the reaction-wise and energy-wise array of RRs/responses
        calculate reaction-wise sum
        """
        spectrum_model_arr = [self.spectrum_model(theta,i) for i in self.energy_grid]
        rr_model_arr = np.array([np.inner(spectrum_model_arr, i) for i in response])
        rr_arr = np.concatenate(rr)
        sigma_rr_arr = np.concatenate(sigma_rr)
        return -0.5 * np.sum(np.log(2 * np.pi * sigma_rr_arr ** 2)
                             +(rr_arr - rr_model_arr) ** 2 / sigma_rr_arr ** 2)

    # def log_likelihood_response_function(self,theta, response, sigma_response, rr):
    #     """ input reaction-wise and energy-wise array of RRs/responses
    #     calculate reaction-wise sum
    #     """
    #     spectrum_model_arr = [self.spectrum_model(theta,i) for i in self.energy_grid]
    #     response_model_arr = np.array([i/spectrum_model_arr for i in rr])
    #     response_arr = np.concatenate(response)
    #     sigma_response_arr = np.concatenate(sigma_response)
    #     return -0.5 * np.sum(np.log(2 * np.pi * sigma_response_arr ** 2)
    #                          +(response_arr - response_model_arr) ** 2 / sigma_response_arr ** 2)

    # posterior for emcee to sample
    def _log_posterior(self,theta,rr,sigma_rr,response):
        """ combined distribution for the measurements (likelihoods)
        and the parameterised neutron spectrum (prior).
        checks for non physical prior first 
        """
        log_prior = np.log(self.spectrum_prior(theta))
        if not np.isfinite(log_prior).any():
            return -np.inf
        return (log_prior + self._log_likelihood(theta,rr,sigma_rr,response))

    # run mcmc
    def run_ensemble_mcmc(self,nparam,param_names,rm_samples,
                          nwalkers,nburn,nsteps,guesses):
        """ run ensemble mcmc
        need to seperate this out into lots of user-controllable functions

        Parameters
        ----------
        isotope_activity : list[float]
            Activity and relative uncertainty of the isotope in Bqe

        Returns
        -------
        correction_factor : float
            correction factor to scale your Faraday cup current by
            to get current on the lithium target
        """
        # filter warnings and random number generator
        warnings.filterwarnings("ignore")
        np.random.seed(42)

        # example data for unfolding
        rr = self.reaction_rates
        rr_u = self.reaction_rate_uncerts
        rm = self.response_matrix
        rm_u = self.response_matrix_uncerts

        # get a gaussian-distributed set of response matrices 
        rm_u_norm = [i*j for i,j in zip(rm,rm_u)]
        rm_dist = np.array([])
        rm_dist = [([np.random.normal(i,j) 
                             for i,j in zip(rm,rm_u_norm)]) 
                             for k in range(rm_samples)]

        # stack the initial parameter guesses as 
        # starting positions for each walker 
        new_guesses = []
        for i in range(nparam):
            new_guesses.append(np.random.randint(low=guesses[i][0]*1e2,
                                                 high=guesses[i][1]*1e2,
                                                 size=nwalkers)/1e2)
        new_guesses = np.stack(new_guesses,axis=1)

        # call the sampler for each response matrices sample 
        # response matrices will vary within their uncertainties
        # the reaction rate distribution is the likelihood
        for i in range(rm_samples):
            print(f'running sampler for response matrices sample {i+1}')
            sampler = emcee.EnsembleSampler(nwalkers, nparam, self._log_posterior,
                                            args=[rr,rr_u,rm_dist[i]])

            # run the sampler for nsteps
            sampler.run_mcmc(new_guesses,nsteps)

        # get parameter results, removing nburn steps 
        samples = sampler.get_chain(flat=True,discard=nburn)
        for i in range(nparam):
            param_ave = (float(np.mean(samples[:,i])))
            param_std = (float(np.std(samples[:,i])))
            print(f'{param_names[i]} = {param_ave:.3f} +- {param_std:.3f}')

        # do corner plotting
        fig = corner.corner(samples,
                            bins=40,
                            labels=param_names)
        plt.savefig('corner')

        # print acceptance fraction and autocorr time
        print("Mean acceptance fraction: {0:.3f}".format(
              np.mean(sampler.acceptance_fraction)))
        print("Mean autocorrelation time: {0:.3f} steps".format(
              np.mean(sampler.get_autocorr_time())))
