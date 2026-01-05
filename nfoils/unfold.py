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
    """ class to do mcmc unfolding with the emcee package
    """

    def __init__(self,files_json_path): 
        """ Initialise Bayesian unfolding class 
        set all unfolding data as numpy arrays of floats

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

        # set group structure attribute and convert to MeV
        # and remove lowest energy value so RF/flux values line up
        group_structure_path = files_paths["group_structure"]
        self.group_structure = (np.fromfile(group_structure_path, sep=","))
        self.group_structure = self.group_structure[1:] * 1e-6

        # set reaction rate attributes
        with open(files_paths["reaction_rates"]) as file:
            reaction_rate_csv = csv.reader(file,delimiter=',')
            self.reaction_rates = np.array(list(reaction_rate_csv),
                                           dtype=float)
        with open(files_paths["reaction_rate_uncerts"]) as file:
            reaction_uncerts_csv = csv.reader(file,delimiter=',')
            self.reaction_rate_uncerts = np.array(list(reaction_uncerts_csv),
                                                  dtype=float)

        # set response matrix attributes
        with open(files_paths["response_matrix"]) as file:
            response_matrix_csv = csv.reader(file,delimiter=',')
            self.response_matrix = np.array(list(response_matrix_csv),
                                            dtype=float)
        with open(files_paths["response_matrix_uncerts"]) as file:
            response_uncerts_csv = csv.reader(file,delimiter=',')
            self.response_matrix_uncerts = np.array(list(response_uncerts_csv),
                                                    dtype=float)

        # FOR TESTING:
        # adjust to unfold just specific parts of the p-li spectrum
        # use [162:] gs values and [6:] rr/rf for just 14 MeV peak
        # use [140:] gs values and [3:] rr/rf for 14 and 8 MeV peak
        # use [70:] gs vals and [2:] rr/rf for 0.1MeV+
        #self.group_structure = self.group_structure[70:]
        #self.reaction_rates = self.reaction_rates[2:]
        #self.reaction_rate_uncerts = self.reaction_rate_uncerts[2:]
        #self.response_matrix = [i[70:] for i in self.response_matrix[2:]]
        #self.response_matrix_uncerts = [i[70:] for i in 
        #                                self.response_matrix_uncerts[2:]]

    def gaussian(self,mean,sigma,peak,energy):
        """ gaussian distribution for defining the model in MeV
        
        Parameters
        ----------
        mean : float
            mean of the distribution
        sigma : float
            width of the distribution
        peak : float
            peak of the distribution
        energy : float
            neutron energy

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        diff = np.sum(energy - mean)
        non_scaled_peak = ( 1 / np.sqrt(2 * np.pi * sigma ** 2))
        scale = 1e-2 * peak / non_scaled_peak
        flux = ((scale / np.sqrt(2 * np.pi * sigma ** 2))
                * np.exp( - 0.5 * diff ** 2 / sigma ** 2))
        return flux

    def lognormal(self,mean,sigma,peak,energy):
        """ lognormal distribution for defining the model in MeV
        
        Parameters
        ----------
        mean : float
            mean of the distribution
        sigma : float
            width of the distribution
        peak : float
            peak of the distribution
        energy : float
            neutron energy

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        lndiff = np.sum(np.log(energy) - mean)
        front_term = np.sum(np.dot(energy, np.sqrt(2 * np.pi * sigma ** 2)))
        non_scaled_peak = ( 1 / front_term)
        scale = 1e-2 * peak / non_scaled_peak
        flux = ((scale / front_term)
                * np.exp( - 0.5 * lndiff ** 2 / sigma ** 2))
        return flux

    def model(self,theta,energy):
        """ model for the neutron flux, given energy and parameters theta. 
        used to create the neutron spectrum prior. 
        needs filling in by the user using the same parameters/returns

        Parameters
        ----------
        theta : tuple
            parameters
        energy : float
            neutron energy

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        print('please create flux model first')
        exit()

    def prior(self,theta):
        """ prior for the neutron spectrum, given parameters theta. 
        combined with the likelihood to create the posterior. 
        needs filling in by the user using the same
        parameters/returns

        Parameters
        ----------
        theta : tuple
            parameters

        Returns
        -------
        prior : float
            prior distribution for the neutron spectrum
            (neutron flux model for the entire group structure)
        """
        print('please create flux prior first')
        exit()

    def _log_likelihood(self,theta, rr, sigma_rr, response):
        """ gaussian distribution for the reaction rate measurements, 
        combined with the prior to create the posterior

        Parameters
        ----------
        theta : tuple
            parameters
        rr : array[array]
            array of reaction rates
        sigma_rr : array[array]
            array of reaction rate uncertainties
        response : array[array]
            array of response functions for each reaction

        Returns
        -------
        log_likelihood : float
            log likelihood for reaction rate measurements
        """
        flux_model_arr = [self.model(theta,i) for i in self.group_structure]
        rr_model_arr = np.array([np.inner(flux_model_arr, i) 
                                 for i in response])
        rr_arr = np.concatenate(rr)
        sigma_rr_arr = np.concatenate(sigma_rr)
        likelihood = -0.5 * np.sum(np.log(2 * np.pi * sigma_rr_arr ** 2)
                     +(rr_arr - rr_model_arr) ** 2 / sigma_rr_arr ** 2)
        return likelihood

    # try and get this to work to avoid doing expensive MC 
    # over the MCMC process to incorporate response function uncertainty
    # def log_likelihood_response_function(self,theta, response, 
    #                                      sigma_response, rr):
    #     """ input reaction-wise and energy-wise array of RRs/responses
    #     calculate reaction-wise sum
    #     """
    #     spectrum_model_arr = [self.spectrum_model(theta,i) 
    #                           for i in self.grop_structure]
    #     response_model_arr = np.array([i/spectrum_model_arr for i in rr])
    #     response_arr = np.concatenate(response)
    #     sigma_response_arr = np.concatenate(sigma_response)
    #     return (-0.5 * np.sum(np.log(2 * np.pi * sigma_response_arr ** 2)
    #             +(response_arr - response_model_arr) ** 2
    #             / sigma_response_arr ** 2))

    def _log_posterior(self,theta,rr,sigma_rr,response):
        """ combined distribution for the measurements (likelihoods)
        and the parameterised neutron spectrum (prior), for emcee to sample
        checks for non physical prior first

        Parameters
        ----------
        theta : tuple
            parameters
        rr : array[array]
            array of reaction rates
        sigma_rr : array[array]
            array of reaction rate uncertainties
        response : array[array]
            array of response functions for each reaction

        Returns
        -------
        log_posterior : float
            log posterior for neutron flux model, given
            the reaction rate measurement disitribution
        """
        log_prior = np.log(self.prior(theta))
        if not np.isfinite(log_prior).any():
            return -np.inf
        log_posterior = (log_prior + 
                         self._log_likelihood(theta,rr,sigma_rr,response))
        return log_posterior

    # maybe seperate this out into some more functions?
    def run_ensemble_mcmc(self,nparam,param_names,rm_samples,
                          nwalkers,nburn,nsteps,guesses):
        """ run ensemble mcmc for given random samples taken from the 
        response matrix distributions (so, MCMCMC). 
        produces corner plot for the parameters - use this to analyse 
        the quality of your model/prior/starting-guesses, and 
        then try again with different inputs. 
        should be an iterative process

        Parameters
        ----------
        nparam : int
            number of parameters in the neutron flux model
        param_names : list[str]
            list of the names of the parameters
        rm_samples : int
            number of samples to take from the response matrix
            distributions. start with 1
        nwalkers : int
            number of MCMC walkers/chains
        nburn : int
            "burn-in" period to let chains stabilize
        nsteps : int
            total number of MCMC steps to take (including nburn)
            no. of trace results =  nwalkers * (nsteps-nburn)
        guesses : list[int]
            set initial guesses for each parameter
            sets starting walker positions

        Returns
        -------
        param_aves : list[float]
            mean values for the parameters
        param_stds : list[float]
            standard deviation values for the parameters
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
            low_guess = guesses[i] - 0.1*guesses[i]
            high_guess = guesses[i] + 0.1*guesses[i]
            new_guesses.append(np.random.randint(low=low_guess*1e2,
                                                 high=high_guess*1e2,
                                                 size=nwalkers)/1e2)
        new_guesses = np.stack(new_guesses,axis=1)

        # call the sampler for each response matrices sample 
        # response matrices will vary within their uncertainties
        # the reaction rate distribution is the likelihood
        for i in range(rm_samples):
            print(f'running sampler for response matrices sample {i+1}')
            sampler = emcee.EnsembleSampler(nwalkers, nparam,
                                            self._log_posterior,
                                            args=[rr,rr_u,rm_dist[i]])

            # run the sampler for nsteps
            sampler.run_mcmc(new_guesses,nsteps)

        # get parameter results, removing nburn steps 
        samples = sampler.get_chain(flat=True,discard=nburn)
        param_aves = []
        param_stds = []
        for i in range(nparam):
            param_ave = (float(np.mean(samples[:,i])))
            param_aves.append(param_ave)
            param_std = (float(np.std(samples[:,i])))
            param_stds.append(param_std)
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
        
        # return parameter values
        return param_aves,param_stds
    
    def plot_spectrum(self,param_guesses,param_aves,param_stds):
        """ plot spectrum with uncertainty using the parameter mean values 
        and the parameter standard deviations

        Parameters
        ----------
        param_guesses : list[float]
            initial guesses for the parameters fed into mcmc
        param_aves : list[float]
            mean values for the parameters
        param_stds : list[float]
            standard deviation values for the parameters
        """
        # convert lists -> arrays and calculate max/min parameters
        param_ave_array = np.array(param_aves)
        param_std_array = np.array(param_stds)
        param_max_array = param_ave_array + param_std_array
        param_min_array = param_ave_array - param_std_array

        # convert param arrays to tuples
        theta = tuple(param_ave_array)
        theta_max = tuple(param_max_array)
        theta_min = tuple(param_min_array)
        theta_guesses = tuple(param_guesses)

        # calculate mean, max and min spectra
        # and the spectrum from the initial guess parameters
        mean_spectrum = [self.model(theta,i) 
                         for i in self.group_structure]
        max_spectrum = [self.model(theta_max,i) 
                        for i in self.group_structure]
        min_spectrum = [self.model(theta_min,i) 
                        for i in self.group_structure]
        guess_spectrum = [self.model(theta_guesses,i) 
                          for i in self.group_structure]

        # plot solution and prior spectrum guess
        fig, ax = plt.subplots(1,figsize=(8,6))
        ax.step(self.group_structure, guess_spectrum,
                label='Initial guess',c='blue',where='pre')
        ax.step(self.group_structure, mean_spectrum,
                label='Solution',c='magenta',where='pre')
        ax.fill_between(self.group_structure, min_spectrum,max_spectrum,
                        alpha=0.25,step='pre')
        #ax.set_xscale('log')
        ax.set_xlim(1,20)
        ax.set_xlabel('Neutron energy (MeV)')
        ax.set_yscale('log')
        ax.set_ylim(1e3,1e8)
        ax.set_ylabel('Flux (n cm$^{-2}$ s$^{-1}$)')
        ax.grid()
        ax.legend(loc="upper left",frameon=True, fontsize=18,fancybox=False,facecolor='white')
        plt.savefig('spectrum.png')
