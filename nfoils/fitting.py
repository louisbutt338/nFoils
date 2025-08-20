""" module for doing fitting of gamma spec efficiency curves
"""

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},weight='normal',size=20)
from scipy.optimize import curve_fit 
import json

class CurveFitter:
    """ class for fitting efficiency functions 
    """
    def __init__(self,input_data_path,input_data_filename,
                interpolation_range_start,interpolation_range_end,
                no_of_monte_carlo_samples):
        """ initialise class
        """

        # set parameters
        self.input_data_path = input_data_path
        self.input_data_filename = input_data_filename
        self.int_range_start = interpolation_range_start
        self.int_range_end = interpolation_range_end
        self.no_of_monte_carlo_samples = no_of_monte_carlo_samples
        
        # change some other parameters
        self.experimental_data = json.load(
            open(f'{input_data_path}/{input_data_filename}.json'))
        self.x_data = [float(i) for i in self.experimental_data.keys()]
        self.y_data = [self.experimental_data[i]["efficiency" ] 
                       for i in self.experimental_data.keys()]
        self.errors = [self.experimental_data[i]["uncertainty"] 
                       for i in self.experimental_data.keys()]
        self.interpolation_range = np.arange(interpolation_range_start,
                                             interpolation_range_end,1)

    def _spec_function(self,energy,a0,a1,a2,a3):
        """ define efficiency polynomial function - output float list

        Parameters
        ----------
        energy : float list
            Gamma energies in keV
        a0 : float
            Zero term parameter
        a1 : float
            1st term parameter for log(E)^1
        a2 : float
            1st term parameter for log(E)^2
        a3 : float
            1st term parameter for log(E)^3
        """
        polynomial = (a0 + a1*np.log(energy)**1 
                      + a2*np.log(energy)**2 
                      + a3*np.log(energy)**3 )
        poly_array = np.exp(polynomial)
        return poly_array

    def _residuals(self,x_list,params_list,y_list):
        """calculated fitted data and find the residuals and rChi2
        """
        fit_data = np.array(self._spec_function(x_list, *params_list))
        residuals = y_list - fit_data 
        chi_squared = np.sum((residuals / self.errors) ** 2)
        dof = len(y_list) - len(params_list) 
        reduced_chi_squared = chi_squared / dof
        return residuals,reduced_chi_squared

    def _single_fit(self,x_list,y_list,param_list):
        """ fit the data once and return the equation params
        """
        params, covs  = curve_fit(self._spec_function, 
                                  x_list, y_list, 
                                  param_list,sigma=self.errors,
                                  absolute_sigma=True)
        errs = np.sqrt(np.diag(covs))
        residuals,reduced_chi_squared = self._residuals(x_list,params,y_list)
        return params,errs,residuals,reduced_chi_squared

    def _monte_carlo_fit(self):
        """ fit the data with MC method
        return the mean, stdev and residuals of the final fit
        """
        N = self.no_of_monte_carlo_samples
        a_samples = []
        a1_samples = []
        a2_samples = []
        a3_samples = []
        mc_solutions = []

        for i in range(N):
            y_mc = self.y_data +np.random.normal(size=len(self.y_data),
                                                 scale=self.errors )
            init_params = self._single_fit(self.x_data,
                                           self.y_data,
                                           [0,0,0,0])[0]
            single_params = self._single_fit(self.x_data,
                                             y_mc,
                                             init_params)[0]
            a0_mc, a1_mc,a2_mc,a3_mc = single_params

            #calculated fitted data and find a solution to average over
            fit_data_mc = self._spec_function(self.interpolation_range,
                                              *single_params)
            a_samples.append (a0_mc)
            a1_samples.append(a1_mc)
            a2_samples.append(a2_mc)
            a3_samples.append(a3_mc)
            mc_solutions.append(fit_data_mc)
              
        # Compute mean parameter values and uncertainties
        a_mc =  np.mean(a_samples)
        a1_mc = np.mean(a1_samples)
        a2_mc = np.mean(a2_samples)
        a3_mc = np.mean(a3_samples)
        params_mc = a_mc,a1_mc,a2_mc,a3_mc

        #calculated fitted data for montecarlo and find the chi squared for MC
        residuals,r_chi_squared = self._residuals(self.x_data,
                                                  params_mc,
                                                  self.y_data)
        
        #return mc_solutions_mean,mc_solutions_std_dev,residuals
        return params_mc,residuals,r_chi_squared,mc_solutions

    def _mc_plotter(self,solutions,standard_dev,residuals):
        """ plot the monte carlo data

        Parameters
        ----------
        fitting_results: arrays,probably
            The output of the monte_carlo_fit function
        """
        # plot function with eror bar
        plt.scatter(self.x_data, self.y_data,
                    label="Data", c='r',marker='o',lw=2)
        plt.plot(self.interpolation_range, solutions,
                 label="Fitted Curve", color='blue')
        plt.fill_between(self.interpolation_range,
                         solutions-standard_dev, solutions+standard_dev,
                         step='post', alpha=0.25)
        plt.errorbar(self.x_data, self.y_data, yerr=self.errors,
                     lw=2,capsize=2,color='k',zorder=-1,fmt='none')
        plt.legend()
        plt.xlim(0,2000)
        plt.ylim(-0.05,0.3)
        plt.xlabel("Gamma energy (keV)")
        plt.ylabel("Efficiency")
        plt.savefig("mc_function.png")
        plt.close()

        #plot MC residuals 
        plt.scatter( self.x_data, residuals,c='r',marker='o',lw=2)
        plt.errorbar(self.x_data, residuals,yerr=self.errors,lw=2,
                     capsize=2,color='k',zorder=-1,fmt='none')
        plt.title("Plot of the residual of the fit")
        plt.xlabel("Time (s)")
        plt.ylabel("Residual")
        plt.savefig("mc_residuals.png")
        plt.close()

    def run_single(self):
        """ run the single fit only
        """
        single_fit_results = self._single_fit()
        a0, a1,a2,a3 = single_fit_results[0]
        a0_err, a1_err,a2_err,a3_err = single_fit_results[1]
        reduced_chi_squared = single_fit_results[2]
        print(f"Estimated Single Fit Parameters: \n a0 = {a0}+/-{a0_err},"
              f" a1 = {a1}+/-{a1_err}, a2 = {a2}+/-{a2_err}, a3 = {a3}+/-{a3_err}"
              f" \n rChi2 = {reduced_chi_squared}")

    def run_mc(self):
        """ run the mc fit
        """
        params_mc,residuals,r_chi_squared,mc_solutions = self._monte_carlo_fit()
        a_mc,a1_mc,a2_mc,a3_mc = params_mc
        print(f"Estimated MC Parameters: \n a0 = {a_mc}, a1 = {a1_mc} ",
              f"a2 = {a2_mc}, a3 = {a3_mc} \n rChi2 = {r_chi_squared} ")

        # calculate the average fit for plotting and the error
        mc_solutions_mean = np.mean(mc_solutions,axis=0)
        mc_solutions_std_dev = np.std(mc_solutions,axis=0)
        mc_frac_uncert = mc_solutions_std_dev/mc_solutions_mean
        print('fractional uncertainty along interpolation range at specified energy '
              f'is {np.mean(mc_frac_uncert[self.int_range_start:self.int_range_end])}')
        
        #plot
        self._mc_plotter(mc_solutions_mean,mc_solutions_std_dev,residuals)