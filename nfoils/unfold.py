""" 
module for doing Bayesian unfolding with emcee

D Foreman-Mackey et al, emcee: The MCMC Hammer, 
Instrumentation and Methods for Astrophysics 2012
"""

import emcee
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=20)

class BayesianUnfolding:
    """ class to do mcmc unfolding 
    """

    def __init__(self,group_structure): 
        """ Initialise TargetAnalysis class

        Attributes
        ----------
        target_json_name : str
            name of the target data json
        """

        # set attributes
        self.energy_grid = group_structure

        # load the target data


    # distributions in straight line example:
    def log_prior_line(self,theta):
        alpha, beta, sigma = theta
        if sigma < 0:
            return -np.inf  # log(0)
        else:
            return -1.5 * np.log(1 + beta ** 2) - np.log(sigma)
        
    def log_likelihood_line(self,theta, x, y):
        alpha, beta, sigma = theta
        y_model = alpha + beta * x
        return -0.5 * np.sum(np.log(2 * np.pi * sigma ** 2) + (y - y_model) ** 2 / sigma ** 2)

    def log_posterior_line(self,theta,x,y):
        return self.log_prior(theta) + self.log_likelihood(theta, x, y)


    # distributions in unfolding example:
    def prior_spectrum(self,theta): # should only be theta here
        mean,sigma,peak = theta
        return -0.5 * peak (2 * np.pi * sigma ** 2) * np.exp((self.energy_grid - mean) ** 2 / sigma ** 2)

    # include sum over energy-bins and reaction-channels in both likelihoods
    def log_likelihood_reaction_rate(self,theta, rr, sigma_rr, response):
        rr_model = self.prior_spectrum(theta) * response
        return -0.5 * np.sum(np.log(2 * np.pi * sigma_rr ** 2) + (rr - rr_model) ** 2 / sigma_rr ** 2)

    def log_likelihood_response_function(self,theta, response, sigma_response, rr):
        response_model = rr / self.prior_spectrum(theta)
        return -0.5 * np.sum(np.log(2 * np.pi * sigma_response ** 2) + (response - response_model) ** 2 / sigma_response ** 2)

    # start simple with gaussian, get more complex later
    def log_posterior_spectrum(self,theta,energy,rr,sigma_rr,response,sigma_response):
        return (np.log(self.prior_spectrum(theta,energy)) 
                + self.log_likelihood_reaction_rate(theta,rr,sigma_rr,response)
                + self.log_likelihood_response_function(theta,response,sigma_response,rr))


    # distribution for multivariate gaussian example
    # searching for theta, given we have some input mu and covs
    def log_posterior(self,theta,mu,cov):
        diff = theta - mu
        return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))


    # run mcmc
    def run_ensemble_mcmc(self):
        """ run the thing for be7 in a lithium target and print results

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

        # try for 5 dimensions (so 5 parameters in the model)
        ndim = 5

        # set up example true data for means (mu) and cov (cov)
        np.random.seed(42)
        means = np.random.rand(ndim)
        print(f'mu values: {means}')
        cov = 0.5 - np.random.rand(ndim**2).reshape((ndim, ndim))
        cov = np.triu(cov)
        cov += cov.T - np.diag(cov.diagonal())
        cov = np.dot(cov, cov)

        # Here we'll set up the computation. emcee combines multiple "walkers",
        # each of which is its own MCMC chain. The number of trace results will
        # be nwalkers * nsteps
        nwalkers = 50  # number of MCMC walkers
        nburn = 100  # "burn-in" period to let chains stabilize
        nsteps = 5000  # number of MCMC steps to take in full simulation

        # set initial guesses for the walker positions in each dimension
        # and call the sampler
        starting_guesses = np.random.rand(nwalkers, ndim)
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_posterior
                                        ,args=[means,cov])

        # run the sampler for burn in
        print('burn-in run')
        state = sampler.run_mcmc(starting_guesses, nburn)
        sampler.reset()

        # run the sampler for production
        print('production run')
        sampler.run_mcmc(state,nsteps)

        # plotting and get aves
        samples = sampler.get_chain(flat=True)
        param_aves = []
        for i in np.arange(ndim):
            param_aves.append(float(np.mean(samples[:,i])))
            plt.figure(figsize=(10,8))
            plt.hist(samples[:, i], 100, color="k", histtype="step")
            plt.xlabel(r"$\theta$")
            plt.ylabel(r"$p(\theta)$")
            plt.gca().set_yticks([])
            plt.savefig(f'param_{i}_dist')
        print(f'theta estimates: {param_aves}')

        # print acceptance fraction and autocorr time
        print("Mean acceptance fraction: {0:.3f}".format(
              np.mean(sampler.acceptance_fraction)))
        print("Mean autocorrelation time: {0:.3f} steps".format(
              np.mean(sampler.get_autocorr_time())))
