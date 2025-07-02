import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},weight='normal',size=20)
from scipy.optimize import curve_fit 
import json

class CurveFitter:
    def __init__(self,input_data_path,input_data_filename,
                interpolation_range_start,interpolation_range_end,
                no_of_monte_carlo_samples):
        self.input_data_path = input_data_path
        self.input_data_filename = input_data_filename
        self.interpolation_range_start = interpolation_range_start
        self.interpolation_range_end = interpolation_range_end
        self.no_of_monte_carlo_samples = no_of_monte_carlo_samples
        # set variables
        self.experimental_data = json.load(open(f'{input_data_path}/{input_data_filename}.json'))
        self.x_data = [float(i) for i in self.experimental_data.keys()]
        self.y_data = [self.experimental_data[i]["efficiency" ] for i in self.experimental_data.keys()]
        self.errors = [self.experimental_data[i]["uncertainty"] for i in self.experimental_data.keys()]
        self.interpolation_range = np.arange(interpolation_range_start,interpolation_range_end,1)

    #define efficiency polynomial
    def _spec_function(self,energy,a0,a1,a2,a3):
        polynomial = a0 + a1*np.log(energy)**1 + a2*np.log(energy)**2 + a3*np.log(energy)**3 
        return np.exp(polynomial)

    #SINGLE FITTING OF DATA
    def _single_fit(self):
        params, covs  = curve_fit(self._spec_function, self.x_data, self.y_data, p0=[0,0,0,0],sigma=self.errors,absolute_sigma=True)
        a0, a1,a2,a3 = params #need this bit to unpack the tuple 
        errs = np.sqrt(np.diag(covs))
        a0_err,a1_err,a2_err,a3_err = errs
        #calculated fitted data and find the residuals etc
        fit_data = self._spec_function(self.x_data, *params)
        residuals = self.y_data - fit_data 
        chi_squared = np.sum((residuals / self.errors) ** 2)
        dof = len(self.y_data) - len(params)  # Degrees of freedom
        reduced_chi_squared = chi_squared / dof
        #print(f"Estimated Single Fit Parameters: \n a0 = {a0}+/-{a0_err}, a1 = {a1}+/-{a1_err}, a2 = {a2}+/-{a2_err}, a3 = {a3}+/-{a3_err} \n rChi2 = {reduced_chi_squared}")
        return params

    #MONTE CARLO FITTING OF DATA
    def _monte_carlo_fit(self):
        N = self.no_of_monte_carlo_samples
        a_samples = []
        a1_samples = []
        a2_samples = []
        a3_samples = []
        a_errs_mc = []
        mc_solutions = [(self._spec_function(self.interpolation_range, *self._single_fit()))]

        for i in range(N):
            y_mc = self.y_data +np.random.normal(size=len(self.y_data), scale=self.errors )
            try:
                params_mc, covs_mc = curve_fit(self._spec_function, self.x_data, y_mc, 
                                               p0=[self._single_fit()[0],self._single_fit()[1],self._single_fit()[2],self._single_fit()[3]],
                                               sigma=self.errors, absolute_sigma=True)
                a0_mc, a1_mc,a2_mc,a3_mc = params_mc #need this bit to unpack the tuple 
                errs_mc = np.sqrt(np.diag(covs_mc))
                a0_err_mc,a1_err_mc,a2_err_mc,a3_err_mc = errs_mc
                #calculated fitted data and find the residuals etc
                fit_data_mc = self._spec_function(self.interpolation_range, *params_mc)
                a_samples.append( a0_mc)
                a1_samples.append(a1_mc)
                a2_samples.append(a2_mc)
                a3_samples.append(a3_mc)
                a_errs_mc.append(errs_mc)
                mc_solutions.append(fit_data_mc)
            except RuntimeError:
                pass  # Skip failed fits
              
        # Compute mean parameter values and uncertainties
        a_mc =  np.mean(a_samples)
        a1_mc = np.mean(a1_samples)
        a2_mc = np.mean(a2_samples)
        a3_mc = np.mean(a3_samples)

        # calculate the average fit for plotting and the error
        mc_solutions_mean = np.mean(mc_solutions,axis=0)
        mc_solutions_std_dev = np.std(mc_solutions,axis=0)
        mc_fractional_uncert = mc_solutions_std_dev/mc_solutions_mean
        print(f'fractional uncertainty along interpolation range at specified energy is {np.mean(mc_fractional_uncert[self.interpolation_range_start:self.interpolation_range_end])}')

        #calculated fitted data for montecarlo and find the chi squared for MC
        fit_data_mc = self._spec_function(self.x_data,  a_mc,a1_mc,a2_mc,a3_mc)
        residuals_mc = self.y_data - fit_data_mc 
        chi_squared_mc = np.sum((residuals_mc / self.errors) ** 2)
        reduced_chi_squared_mc = chi_squared_mc / (len(self.y_data) - len(self._single_fit()))
        print(f"Estimated MC Parameters: \n a0 = {a_mc}, a1 = {a1_mc}, a2 = {a2_mc}, a3 = {a3_mc} \n rChi2 = {reduced_chi_squared_mc} ")
        
        return mc_solutions_mean,mc_solutions_std_dev,residuals_mc

    #plot data
    def _plotter(self,fitting_results):
        solutions = fitting_results[0]
        standard_dev = fitting_results[1]
        residuals = fitting_results[2]
        plt.scatter(self.x_data, self.y_data, label="Data", c='r',marker='o',lw=2)
        plt.plot(self.interpolation_range, solutions, label="Fitted Curve", color='blue')
        plt.fill_between(self.interpolation_range, solutions-standard_dev, solutions+standard_dev, step='post', alpha=0.25)
        plt.errorbar(self.x_data, self.y_data, yerr=self.errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
        plt.legend()
        plt.xlim(0,2000)
        plt.ylim(-0.05,0.3)
        plt.xlabel("Gamma energy (keV)")
        plt.ylabel("Efficiency")
        #plt.savefig("mc_function.png")
        plt.close()

        #plot MC residuals 
        plt.scatter( self.x_data, residuals,c='r',marker='o',lw=2)
        plt.errorbar(self.x_data, residuals,yerr=self.errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
        plt.title("Plot of the residual of the fit")
        plt.xlabel("Time (s)")
        plt.ylabel("Residual")
        #plt.savefig("mc_residuals.png")
        plt.close()

        ## plot distribution of a0 params
        #plt.hist(a_samples, 100)
        #plt.title("Plot of parameter distribution")
        #plt.xlabel("Parameter value")
        #plt.ylabel("Frequency")
        ##plt.savefig("mc_param_dist.png")
        #plt.close()

    def run(self):
        self._plotter(self._monte_carlo_fit())
