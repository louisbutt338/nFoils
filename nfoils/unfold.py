""" 
module for doing Bayesian unfolding with emcee

D Foreman-Mackey et al, emcee: The MCMC Hammer, 
Instrumentation and Methods for Astrophysics 2012
"""

import emcee
import corner
import numpy as np
import csv
import json
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=20)
from scipy import stats

class BayesianUnfolding:
    """ class to do mcmc unfolding with the emcee package
    """

    def __init__(self,files_json_path,nparam,param_names,guesses): 
        """ Initialise Bayesian unfolding class 
        set all unfolding data as numpy arrays of floats

        Attributes
        ----------
        files_json_path : str
            path to the files data json
        nparam : int
            number of parameters in the neutron flux model
        param_names : list[str]
            list of the names of the parameters
        guesses : list[float]
            initial guesses for the parameters fed into mcmc
        """
        # set parameter attributes
        self.nparam = nparam
        self.param_names = param_names
        self.guesses = guesses

        # load the datafiles paths
        with open(files_json_path
                  ) as files_json:
            files_paths = json.load(files_json)

        # set group structure attribute and convert to MeV
        # and remove lowest energy value so RF/flux values line up
        group_structure_path = files_paths["group_structure"]
        self._group_structure_raw = (np.fromfile(group_structure_path, sep=","))
        self._group_structure = self._group_structure_raw[1:] * 1e-6

        # set reaction rate attributes
        with open(files_paths["reaction_rates"]) as file:
            reaction_rate_csv = csv.reader(file,delimiter=',')
            self._reaction_rates = np.array(list(reaction_rate_csv),
                                           dtype=float)
        with open(files_paths["reaction_rate_uncerts"]) as file:
            reaction_uncerts_csv = csv.reader(file,delimiter=',')
            self._reaction_rate_uncerts = np.array(list(reaction_uncerts_csv),
                                                  dtype=float)

        # set response matrix attributes
        with open(files_paths["response_matrix"]) as file:
            response_matrix_csv = csv.reader(file,delimiter=',')
            self._response_matrix = np.array(list(response_matrix_csv),
                                            dtype=float)
        with open(files_paths["response_matrix_uncerts"]) as file:
            response_uncerts_csv = csv.reader(file,delimiter=',')
            self._response_matrix_uncerts = np.array(list(response_uncerts_csv),
                                                    dtype=float)

    # getters and setters for the data attributes

    @property
    def group_structure(self):
        """ gets group structure attribute
        
        Returns
        -------
        group_structure : array
            group structure in MeV
        """
        return self._group_structure

    @group_structure.setter
    def group_structure(self, group_structure):
        """ sets group structure attribute 
        
        Parameter
        ---------
        group_structure : arr[arr]
            group structure in MeV
        """
        self._group_structure = group_structure

    @property
    def reaction_rates(self):
        """ gets reaction rates attribute
        
        Returns
        -------
        reaction_rates : arr[arr]
            foil reaction rates
        """
        return self._reaction_rates

    @reaction_rates.setter
    def reaction_rates(self, reaction_rates):
        """ sets reaction rates attribute 
        
        Parameter
        ---------
        reaction_rates : arr[arr]
            foil reaction rates
        """
        self._reaction_rates = reaction_rates

    @property
    def reaction_rate_uncerts(self):
        """ gets reaction rate uncertainties attribute
        
        Returns
        -------
        reaction_rate_uncerts : arr[arr]
            foil reaction rate uncertainties
        """
        return self._reaction_rate_uncerts

    @reaction_rate_uncerts.setter
    def reaction_rate_uncerts(self, reaction_rate_uncerts):
        """ sets reaction rate uncertainties attribute 
        
        Parameter
        ---------
        reaction_rate_uncerts : arr[arr]
            foil reaction rate uncertainties
        """
        self._reaction_rate_uncerts = reaction_rate_uncerts

    @property
    def response_matrix(self):
        """ gets response matrix attribute
        
        Returns
        -------
        response_matrix : arr[arr]
            foil reaction response matrix
        """
        return self._response_matrix

    @response_matrix.setter
    def response_matrix(self, response_matrix):
        """ sets response matrix attribute 
        
        Parameter
        ---------
        response_matrix : arr[arr]
            foil reaction response matrix
        """
        self._response_matrix = response_matrix

    @property
    def response_matrix_uncerts(self):
        """ gets response matrix uncertainties attribute
        
        Returns
        -------
        response_matrix_uncerts : arr[arr]
            foil reaction response matrix uncertainties
        """
        return self._response_matrix_uncerts

    @response_matrix_uncerts.setter
    def response_matrix_uncerts(self, response_matrix_uncerts):
        """ sets response matrix uncertainties attribute 
        
        Parameter
        ---------
        response_matrix_uncerts : arr[arr]
            foil reaction response matrix uncerts
        """
        self._response_matrix_uncerts = response_matrix_uncerts

    # example distributions to make your model with

    def gaussian(self,mean,sigma,peak,energy):
        """ gaussian distribution 
        see scipy.stats documentation for further details
        
        Parameters
        ----------
        mean : float
            mean of the distribution (MeV)
        sigma : float
            width of the distribution (MeV)
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        non_scaled_peak = stats.norm.pdf(mean,loc=mean,scale=sigma)
        scale = peak / non_scaled_peak
        flux = scale*stats.norm.pdf(energy,loc=mean,scale=sigma)
        return flux
    
    def skewgaussian(self,mean,sigma,peak,skew,energy):
        """ skewed gaussian distribution 
        see scipy.stats documentation for further details
        
        Parameters
        ----------
        mean : float
            mean of the distribution (MeV)
        sigma : float
            width of the distribution (MeV)
        peak : float
            peak of the distribution (n/cm2/s)
        skew : float
            skew of distribution
        energy : float
            neutron energy MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        non_scaled_peak = stats.skewnorm.pdf(mean,skew,loc=mean,scale=sigma)
        scale = peak / non_scaled_peak
        flux = scale*stats.skewnorm.pdf(energy,skew,loc=mean,scale=sigma)
        return flux

    def lognormal(self,mode,sigma,peak,energy):
        """ lognormal distribution 
        see scipy.stats documentation for further details
        
        Parameters
        ----------
        mode : float
            mode of the distribution (MeV)`
        sigma : float
            width of the logarithm of the distribution
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        mu = np.log(mode) + sigma**2
        non_scaled_peak = stats.lognorm.pdf(mode,s=sigma,scale=np.exp(mu))
        scale = peak / non_scaled_peak
        flux= scale * stats.lognorm.pdf(energy,s=sigma,scale=np.exp(mu))
        return flux
    
    def gamma(self,alpha,beta,peak,energy):
        """ gamma distribution 
        see scipy.stats documentation for further details

        Parameters
        ----------
        alpha : float
            shape parameter for distribution
        beta : float
            scale parameter for the distribution
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        mode = (alpha-1)*beta
        non_scaled_peak = stats.gamma.pdf(mode,a=alpha,scale=beta)
        scale=peak/non_scaled_peak
        flux = scale * stats.gamma.pdf(energy,a=alpha,scale=beta)
        return flux
    
    def weibull(self,alpha,beta,peak,energy):
        """ weibull distribution 
        see scipy.stats documentation for further details

        Parameters
        ----------
        alpha : float
            shape parameter for distribution
        beta : float
            scale parameter for the distribution
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy in MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        mode = beta*((alpha-1)/beta)**(1/alpha)
        non_scaled_peak = stats.weibull_min.pdf(mode,c=alpha,scale=beta)
        scale=peak/non_scaled_peak
        flux = scale * stats.weibull_min.pdf(energy,c=alpha,scale=beta)
        return flux

    def maxwellian(self,temp,peak,energy):
        """ simple maxwellian peak distribution
        see www-ucjf.troja.mff.cuni.cz/cnr11/presentations_dir/kawano.pdf

        Parameters
        ----------
        temp : float
            mode/evaporation temperature (MeV) (t>0)
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy in MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        non_scaled_peak = np.sqrt(temp)* np.exp(-1)
        scale = peak/non_scaled_peak
        flux=(scale*np.sqrt(energy)*
              np.exp(-energy/temp))
        return flux
    
    def watt(self,temp,beta,peak,energy):
        """ Watt fission distribution (modified maxwellian)
        see www-ucjf.troja.mff.cuni.cz/cnr11/presentations_dir/kawano.pdf

        Parameters
        ----------
        temp : float
            mode/evaporation temperature (MeV) (t>0)
        beta : float
            scale parameter for the distribution
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy in MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        non_scaled_peak = (np.sinh(np.sqrt(beta*temp))*
                           np.exp(-1))
        scale = peak/non_scaled_peak
        flux = scale * (np.sinh(np.sqrt(beta*energy))*
                        np.exp(-energy/temp))
        return flux
    
    def weisskopf(self,temp,peak,energy):
        """ weisskopf/evaporation peak distribution (modified maxwellian)
        see www-ucjf.troja.mff.cuni.cz/cnr11/presentations_dir/kawano.pdf

        Parameters
        ----------
        temp : float
            mode/evaporation temperature (MeV) (t>0)
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy in MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        non_scaled_peak = temp * np.exp(-1)
        scale = peak/non_scaled_peak
        flux = (scale*energy*np.exp(-energy/temp))
        return flux

    def epithermal(self,alpha,beta,scale,e_limit,energy):
        """ epithermal straight line distribution 
        see FRUIT paper for details (doi.org/10.1016/j.nima.2007.07.033)

        Parameters
        ----------
        alpha : float
            slope parameter 1 (t>0.5)
        beta : float
            slope parameter 2 (0-1)
        scale : float
            scale of the distribution (n/cm2/s)
        e_limit : float
            lower energy limit ~1e-7 (MeV)
        energy : float
            neutron energy in MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        term1 = 1-np.exp(-(energy/e_limit)**2)
        term2 = energy**(alpha-1)
        term3 = np.exp(-energy/beta)
        flux = scale*term1*term2*term3
        return flux
    
    def thermal(self,mean,peak,energy):
        """ thermal peak distribution 
        see FRUIT paper for details (doi.org/10.1016/j.nima.2007.07.033)

        Parameters
        ----------
        mean : float
            mean/thermal neutron energy ~1.53e-8 (MeV)
        peak : float
            peak of the distribution (n/cm2/s)
        energy : float
            neutron energy in MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        non_scaled_peak = (1/mean)*np.exp(-1)
        scale=peak/non_scaled_peak
        flux=(scale*(energy/mean**2)
              *np.exp(-energy/mean))
        return flux
    
    def powerlaw(self,alpha,scale,energy):
        """ simple power law distribution. 
        Use for the low-energy neutron spectrum
        
        Parameters
        ----------
        alpha: float
            exponent/slope of the distribution 1-2
        scale: float
            scale of the distribution (n/cm2/s)
            (hypothetical flux around 1 MeV)
        energy : float
            neutron energy MeV

        Returns
        -------
        flux : float
            neutron flux for the given energy
        """
        flux = scale / (np.power(energy,(1-alpha)))
        return flux
    
    # pre-set probability distributions for likelihood and posterior

    def log_likelihood(self,theta, model,rr, sigma_rr, response):
        """ gaussian distribution for the reaction rate measurements, 
        combined with the prior to create the posterior

        Parameters
        ----------
        theta : tuple
            parameters
        model : callable
            neutron flux model
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
        model_for_gs = [model(theta,i) for i in self.group_structure]
        rr_model_arr = np.array([np.inner(model_for_gs, i) 
                                 for i in response])
        rr_arr = np.concatenate(rr)
        sigma_rr_arr = np.concatenate(sigma_rr)
        loglikelihood = -0.5 * np.sum(np.log(2 * np.pi * sigma_rr_arr ** 2)
                     +(rr_arr - rr_model_arr) ** 2 / sigma_rr_arr ** 2)
        return loglikelihood

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

    def log_posterior(self,theta,model,log_prior,log_likelihood,
                      rr,sigma_rr,response):
        """ combined distribution for the measurements (likelihoods)
        and the parameterised neutron spectrum (prior), for emcee to sample
        checks for non physical prior first

        Parameters
        ----------
        theta : tuple
            parameters
        model : callable
            neutron flux model function
        log_prior : callable
            log prior function
        log_likelihood : callable
            log likelihood function
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
        if not np.isfinite(log_prior(theta,model)).any():
            return -np.inf
        log_posterior = (log_prior(theta,model) + 
                         log_likelihood(theta,model,rr,sigma_rr,response))
        return log_posterior
    
    # run methods

    def _setup_sampler(self,rm_samples,nwalkers):
        """ set up gaussian dsitributed response matrices to sample from 
        and initial parameter guesses for each walker

        Parameters
        ----------
        rm_samples : int
            number of samples to take from the response matrix
            distributions. start with 1
        nwalkers : int
            number of MCMC walkers/chains

        Returns
        -------
        new_guesses : arr[arr]
            starting guesses for each walker
        rm_dist : arr[arr]
            gaussian-distributed set of response matrices
        """

        # random number generator
        np.random.seed(42)

        # get a gaussian-distributed set of response matrices
        rm = self.response_matrix
        rm_u = self.response_matrix_uncerts
        rm_u_norm = [i*j for i,j in zip(rm,rm_u)]
        rm_dist = np.array([])
        rm_dist = [([np.random.normal(i,j) 
                     for i,j in zip(rm,rm_u_norm)]) 
                     for k in range(rm_samples)]

        # stack the initial parameter guesses as 
        # starting positions for each walker 
        new_guesses = []
        for i in range(self.nparam):
            low_guess  = self.guesses[i] - 0.1*self.guesses[i]
            high_guess = self.guesses[i] + 0.1*self.guesses[i]
            new_guesses.append(np.random.randint(low=low_guess*1e8,
                                                 high=high_guess*1e8,
                                                 size=nwalkers)/1e8)
        new_guesses = np.stack(new_guesses,axis=1)

        return new_guesses,rm_dist

    def run_sampler(self,log_post,model,log_prior,log_likelihood,
                    rm_samples,nwalkers,nburn,nsteps,pool=None):
        """ run ensemble mcmc for given random samples taken from the 
        response matrix distributions. The reaction rate distribution 
        is the likelihood in the model.
        Recommend the following method:
        (1) run sampler iteratively using multiple rm_samples but low
        nwalkers/nsteps, to improve your model/prior/guesses. 
        Use the results/corner plot to improve your prior
        (2) once you have a good prior/guesses, increase 
        nwalkers/nsteps as much as possible. Feel free to bastardise
        this module to run it on a HPC, otherwise will be v slow

        Parameters
        ----------
        log_post : callable
            set to the log_posterior function
        model : callable 
            set to the model function
        log_prior : callable
            set to the log prior function
        log_likelihood : callable
            set to the log likelihood function
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
        pool : callable
            set to multiprocessing.Pool() for parallel, or
            defaults to None

        Returns
        -------
        samples : callable
            output chain from ensemble sampler
        """

        # get the reaction rate data
        rr = self.reaction_rates
        rr_u = self.reaction_rate_uncerts

        # get the guesses for each walker and response matrix distributions
        new_guesses,rm_dist = self._setup_sampler(rm_samples,nwalkers)

        # call the sampler for each response matrix sample 
        # response functions will vary within their uncertainties
        # the reaction rate distribution is the likelihood
        all_samples = []
        for i in range(rm_samples):
            print(f'running sampler for response matrix sample {i+1}')
            sampler = emcee.EnsembleSampler(nwalkers, self.nparam,
                                            log_post, 
                                            args=[model,log_prior,
                                                  log_likelihood,
                                                  rr,rr_u,rm_dist[i]],
                                            pool=pool)
            # run the sampler for nsteps
            sampler.run_mcmc(new_guesses,nsteps)

            # print acceptance fraction and autocorr time
            print("mean acceptance fraction: {0:.3f}".format(
                  np.mean(sampler.acceptance_fraction)))
            #print("mean autocorrelation time: {0:.3f} steps".format(
            #      np.mean(sampler.get_autocorr_time())))

            # get sampler results and append, removing nburn steps 
            samples = sampler.get_chain(flat=True,discard=nburn)
            all_samples.append(samples)

        # join sample results together 
        concat_samples = np.concatenate(all_samples)

        return concat_samples
    
    # postprocessing methods

    def postpro_sampler(self,samples,plotname='corner'):
        """ take results (samples) from run_mcmc. generate 
        a corner plot and output the values for parameters and
        their covariance matrix.
        use corner plot to analyse the quality of your model/prior 
        and then try again with different inputs - iterative process

        Parameters
        ----------
        samples : callable
            output chain from the ensemble sampler
        plotname : str
            name of the corner plot

        Returns
        -------
        param_aves : list[float]
            mean values for the parameters
        param_cov_matrix : list[float]
            covariance matrix for the parameters
        """

        # get parameters - not stds
        param_aves = []
        #param_stds = []
        print('Parameter results:')
        for i in range(self.nparam):
            param_ave = (float(np.mean(samples[:,i])))
            param_aves.append(param_ave)
            #param_std = (float(np.std(samples[:,i])))
            #param_stds.append(param_std)
            print(f'{self.param_names[i]} = {param_ave:.3f}')

        # get cov matrix
        param_cov_matrix = np.cov(samples, rowvar=False)

        # do corner plotting
        fig = corner.corner(samples,
                            bins=40,
                            labels=self.param_names)
        plt.savefig(f'{plotname}.png')
        
        return param_aves,param_cov_matrix
    
    def lethargy_conv(self,spectrum,group_structure):
        """ convert a spectrum into lethargy spectrum

        Parameters
        ----------
        spectrum : list[float]
            spectrum in n/cm2/s (length n)
        group_structure : list[float]
            group structure in MeV (length n+1)
        
        Returns
        -------
        spectrum_leth : list[float]
            the lethargy specturm
        """
        spectrum_leth = []
        for i,j in enumerate(spectrum):
            e_u = group_structure[i+1]
            e_l = group_structure[i]
            leth = np.log(e_u/e_l)
            spectrum_leth.append(j/leth)
        return spectrum_leth
    
    def get_spectra(self,model,param_aves,param_stds,cutoff=None):
        """ get spectra from the model params

        Parameters
        ----------
        model : callable 
            set to the model function
        param_aves : list[float]
            mean values for the parameters
        param_stds : list[float]
            standard deviation values for the parameters
        cutoff : int
            number of unphysical values 
            to cut off from the end of the spectrum

        Returns
        -------
        mean_spectrum : list[float]
            the solution spectrum
        max_spectrum : list[float]
            the upper bound spectrum
        min_spectrum : list[float]
            the lower bound spectrum
        guess_spectrum : list[float]
            the prior guess spectrum
        """
        # shorthand for the unfolding group structure
        gs = self.group_structure

        # calculate max/min parameters and convert all to tuples
        param_ave_array = np.array(param_aves)
        param_std_array = np.array(param_stds)
        theta_max = tuple(param_ave_array + param_std_array)
        theta_min = tuple(param_ave_array - param_std_array)
        theta_guesses = tuple(self.guesses)

        # calculate mean, max, min, initial guess spectra
        max_spectrum = [model(theta_max,i) for i in gs]
        min_spectrum = [model(theta_min,i) for i in gs]
        guess_spectrum = [model(theta_guesses,i) for i in gs]
        mean_spectrum = [np.mean([i,j]) for i,j 
                         in zip(min_spectrum,max_spectrum)]

        # remove preset unphysical values from spectra
        if cutoff != None:
            for i in range(len(gs)-cutoff,len(gs)):
                mean_spectrum[i]=0
                max_spectrum[i]=0
                min_spectrum[i]=0
                guess_spectrum[i]=0

        return (mean_spectrum,max_spectrum,
                min_spectrum,guess_spectrum)
    
    def plot_simple_spectrum(self,model,param_aves,param_stds,
                             cutoff=None,plotname='spectrum'):
        """ plot solution spectrum+uncertainty and initial
        parameter guess spectrum. simple linear plot

        Parameters
        ----------
        model : callable 
            set to the model function
        param_aves : list[float]
            mean values for the parameters
        param_stds : list[float]
            standard deviation values for the parameters
        cutoff : int
            number of unphysical values 
            to cut off from the end of the spectrum
        plotname : str
            name of the spectrum plot
        """
        # get all the spectra and the group structure
        gs = self.group_structure
        (mean_spectrum,max_spectrum,
        min_spectrum,guess_spectrum) = self.get_spectra(model,
                                                        param_aves,
                                                        param_stds,
                                                        cutoff)

        # plot solution and prior spectrum guess
        fig, ax = plt.subplots(1,figsize=(10,6))
        ax.step(gs, guess_spectrum,
                label='Prior guess',c='blue',where='pre')
        ax.step(gs, mean_spectrum,
                label='Solution',c='magenta',where='pre')
        ax.fill_between(gs, min_spectrum,max_spectrum,
                        alpha=0.25,step='pre')
        #ax.set_xscale('log')
        ax.set_xlim(0)
        ax.set_xlabel('Neutron energy (MeV)')
        ax.set_yscale('log')
        ax.set_ylim(1e0)
        ax.set_ylabel('Flux (n cm$^{-2}$ s$^{-1}$)')
        ax.grid()
        ax.legend(loc="lower left",frameon=True, fontsize=18,
                  fancybox=False,facecolor='white')
        fig.tight_layout()
        plt.savefig(f'{plotname}.png')


    def plot_split_spectrum(self,model,param_aves,param_stds,
                            cutoff=None,plotname='spectrum'):
        """ plot solution spectrum+uncertainty and initial
        parameter guess spectrum. split log-linear plot

        Parameters
        ----------
        model : callable 
            set to the model function
        param_aves : list[float]
            mean values for the parameters
        param_stds : list[float]
            standard deviation values for the parameters
        cutoff : int
            number of unphysical values 
            to cut off from the end of the spectrum
        plotname : str
            name of the spectrum plot
        """

        # get all the spectra and the group structure
        gs = self.group_structure
        (mean_spectrum,max_spectrum,
        min_spectrum,guess_spectrum) = self.get_spectra(model,
                                                        param_aves,
                                                        param_stds,
                                                        cutoff)
        
        # cutoff swithces for the C/T graph
        if cutoff!= None:
            gs_cutoff = gs[:-cutoff]
            ones_cutoff = np.ones(len(gs[:-cutoff]))
        else:
            gs_cutoff = gs
            ones_cutoff = np.ones(len(gs))

        fig, ((ax1,ax2),(ax3,ax4)) = (plt.subplots(2,2,figsize=(13,7),
                                    gridspec_kw={'width_ratios': [1, 2],
                                                 'height_ratios': [4, 1]}))
        # plot a priori, solution and uncertainty 
        ax1.step(gs, guess_spectrum, label='Prior guess',where='pre',c='blue')
        ax2.step(gs, guess_spectrum, label='Prior guess',where='pre',c='blue')
        ax1.step(gs, mean_spectrum, label='Solution', where='pre',c='magenta')
        ax2.step(gs, mean_spectrum, label='Solution', where='pre',c='magenta')
        ax1.fill_between(gs, min_spectrum, max_spectrum, step='pre', alpha=0.25)
        ax2.fill_between(gs, min_spectrum, max_spectrum, step='pre', alpha=0.25)

        # plot C/T graph with uncertainty
        soln_ct = np.array(mean_spectrum)/np.array(guess_spectrum)
        min_ct = np.array(min_spectrum)/np.array(guess_spectrum)
        max_ct = np.array(max_spectrum)/np.array(guess_spectrum)
        ax3.step(gs, soln_ct, where='pre',c='magenta')
        ax4.step(gs, soln_ct, where='pre',c='magenta')
        ax3.fill_between(gs, min_ct, max_ct, step='pre',alpha=0.25)
        ax4.fill_between(gs, min_ct, max_ct, step='pre',alpha=0.25)
        ax3.step(gs_cutoff, ones_cutoff, where='post',c='blue')
        ax4.step(gs_cutoff, ones_cutoff, where='post',c='blue')

        # plotting parameters for all 4 graphs 
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylim(1e0)
        ax1.set_xlim(1e-7,1e0)
        ax1.tick_params(axis='x',labelbottom=False)
        ax1.grid()
        ax1.set_ylabel('Flux (n cm$^{-2}$ s$^{-1}$)',y=0.5)
        ax2.set_yscale('log')
        ax2.set_ylim(1e0)
        ax2.set_xlim(1,15)
        ax2.tick_params(axis='y',labelleft=False)
        ax2.tick_params(axis='x',labelbottom=False)
        ax2.grid()
        ax2.legend(loc="upper left", bbox_to_anchor=(0.03, 0.263), borderaxespad=0,
                   frameon=True, fontsize=18,fancybox=False,facecolor='white',framealpha=1)
        ax3.set_xscale('log')
        ax3.set_yscale('log')
        ax3.set_ylim(1e-1,1e1)
        ax3.set_xlim(1e-7,1)
        ax3.grid()
        ax3.set_ylabel('C/M',y=0.5)
        ax4.set_yscale('log')
        ax4.set_ylim(1e-1,1e1)
        ax4.set_xlim(1,15)
        ax4.tick_params(axis='y',labelleft=False)
        ax4.grid()
        fig.supxlabel('Neutron energy (MeV)',y=0.04)
        fig.tight_layout()
        plt.show
        plt.subplots_adjust(wspace=0.04, hspace=0.1)
        plt.savefig(f'{plotname}.png')

    def dump_spectrum(self,model,param_aves,param_stds,
                      cutoff=None,txtname='spectrum'):
        """ dump spectrum for further analysis

        Parameters
        ----------
        model : callable 
            set to the model function
        param_aves : list[float]
            mean values for the parameters
        param_stds : list[float]
            standard deviation values for the parameters
        cutoff : int
            number of unphysical values 
            to cut off from the end of the spectrum
        txtname : str
            name of the dumped spectrum txt file
        """
        # get all the spectra and the group structure
        (mean_spectrum,max_spectrum,
        min_spectrum,guess_spectrum) = self.get_spectra(model,
                                                        param_aves,
                                                        param_stds,
                                                        cutoff)

        # get raw uncertainty and dump spectrum+uncert
        print(f'total flux = {np.sum(mean_spectrum)} n/cm2/s')
        uncert_spectrum=[(i-j) for i,j in zip(max_spectrum,mean_spectrum)]
        np.savetxt(f'{txtname}.txt', 
                   np.transpose([mean_spectrum,uncert_spectrum]))

    def residual_sum_squares(self,spectrum):
        """ find residual sum of sqaures for reaction rates, for 
        a given converged unfolding model

        
        Parameters
        ----------
        spectrum : list[float]
            final unfolded spectrum

        Returns
        -------
        rss : float
            residual sum of squares for unfolding model
        """
        pred_rates = [i*spectrum for i in self.response_matrix]
        obs_rates = self.reaction_rates
        squares_list = [(i-j)**2 for i,j in zip(obs_rates,pred_rates)]
        rss = np.sum(squares_list)
        return rss

    def get_bic(self,model,param_aves,param_stds,cutoff=None):
        """ calculate bayesian information criterion, 
        for a given converged unfolding model

        Parameters
        ----------
        model : callable 
            set to the model function
        param_aves : list[float]
            mean values for the parameters
        param_stds : list[float]
            standard deviation values for the parameters
        cutoff : int
            number of unphysical values 
            to cut off from the end of the spectrum

        Returns
        -------
        bic : float
            Bayesian information criterion for model
        """
        mean_spectrum= self.get_spectra(model,param_aves,
                                        param_stds,cutoff)[0]
        rss = self.residual_sum_squares(mean_spectrum)
        num_rrs = len(self.reaction_rates) 
        bic = (num_rrs*np.log(rss/num_rrs) + self.nparam*np.log(num_rrs))
        return bic