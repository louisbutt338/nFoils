import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
from scipy.optimize import curve_fit # type: ignore
import scipy.optimize as opt # type: ignore
import json

###########################
####### USER INPUTS #######

#experimental data
b03_endcap_exp_data = {
  121.7817 : {"efficiency" : 0.119797036, "uncertainty": 0.006026234 },
  244.6974 : {"efficiency" : 0.091639719, "uncertainty": 0.004623184 },
  344.2785 : {"efficiency" : 0.069565395, "uncertainty": 0.003499949 },
  443.9606 : {"efficiency" : 0.058245795, "uncertainty": 0.002954337 },
  778.9045 : {"efficiency" : 0.035923606, "uncertainty": 0.001811229 },
  867.3800 : {"efficiency" : 0.033051461, "uncertainty": 0.001685693 },
  964.0570 : {"efficiency" : 0.03087022 , "uncertainty": 0.001554952 },
  1085.837 : {"efficiency" : 0.029895735, "uncertainty": 0.001513098 },
  1112.076 : {"efficiency" : 0.027864365, "uncertainty": 0.001404838 },
  1408.013 : {"efficiency" : 0.023700377, "uncertainty": 0.001193829 }
}
b03_38cm_exp_data = {
  121.7817 : {"efficiency" : 0.000972425, "uncertainty": 5.00169E-05 },
  244.6974 : {"efficiency" : 0.00074166 , "uncertainty": 4.19527E-05 },
  344.2785 : {"efficiency" : 0.000563953, "uncertainty": 2.94885E-05 },
  443.9606 : {"efficiency" : 0.000471526, "uncertainty": 3.54344E-05 },
  778.9045 : {"efficiency" : 0.00029065 , "uncertainty": 1.69317E-05 },
  867.3800 : {"efficiency" : 0.000267581, "uncertainty": 2.07387E-05 },
  964.0570 : {"efficiency" : 0.000249694, "uncertainty": 1.4535E-05  },
  1085.837 : {"efficiency" : 0.000241937, "uncertainty": 1.49791E-05 },
  1112.076 : {"efficiency" : 0.00022527 , "uncertainty": 1.33982E-05 },
  1408.013 : {"efficiency" : 0.000191711, "uncertainty": 1.09422E-05 }
}
g11_endcap_exp_data = {
  121.7817 : {"efficiency" : 0.146377303, "uncertainty": 0.007362555 },
  244.6974 : {"efficiency" : 0.07887771 , "uncertainty": 0.003978394 },
  344.2785 : {"efficiency" : 0.048097055, "uncertainty": 0.002419776 },
  443.9606 : {"efficiency" : 0.036811079, "uncertainty": 0.001868924 },
  778.9045 : {"efficiency" : 0.017591791, "uncertainty": 0.000887509 },
  867.3800 : {"efficiency" : 0.015047731, "uncertainty": 0.000772899 },
  964.0570 : {"efficiency" : 0.013543158, "uncertainty": 0.000683006 },
  1085.837 : {"efficiency" : 0.013328323, "uncertainty": 0.00067562  },
  1112.076 : {"efficiency" : 0.011142293, "uncertainty": 0.000562644 },
  1408.013 : {"efficiency" : 0.009452859, "uncertainty": 0.000476738 }
}
g11_10cm_exp_data = {
  121.7817 : {"efficiency" : 0.002770881, "uncertainty": 0.000139351},
  244.6974 : {"efficiency" : 0.001490941, "uncertainty": 0.000185056},
  344.2785 : {"efficiency" : 0.000907987, "uncertainty": 8.83649E-05},
  443.9606 : {"efficiency" : 0.000698606, "uncertainty": 0.000222016},
  778.9045 : {"efficiency" : 0.000332174, "uncertainty": 7.03765E-05},
  867.3800 : {"efficiency" : 0.000278539, "uncertainty": 0.000121242},
  964.0570 : {"efficiency" : 0.000258369, "uncertainty": 5.62885E-05},
  1085.837 : {"efficiency" : 0.000231966, "uncertainty": 7.26431E-05},
  1112.076 : {"efficiency" : 0.000202716, "uncertainty": 5.42226E-05},
  1408.013 : {"efficiency" : 0.000164661, "uncertainty": 3.66736E-05}
}

#select dataset
input_data_path = f"/Users/ljb841@student.bham.ac.uk/nFoils/data/gamma_spec_efficiencies"
input_data_filename = "b03_38cm_data"

# energy range in keV for fitting the data in montecarlo and finding the average uncertainty
interpolation_range_start = 100
interpolation_range_end = 1800

###########################
###########################

# set variables
experimental_data = json.load(open(f'{input_data_path}/{input_data_filename}.json'))
x_data = [float(i) for i in experimental_data.keys()]
y_data = [experimental_data[i]["efficiency" ] for i in experimental_data.keys()]
errors = [experimental_data[i]["uncertainty"] for i in experimental_data.keys()]
interpolation_range = np.arange(interpolation_range_start,interpolation_range_end,1)

#define efficiency polynomial
def spec_function(energy,a0,a1,a2,a3):
    polynomial = a0 + a1*np.log(energy)**1 + a2*np.log(energy)**2 + a3*np.log(energy)**3 
    return np.exp(polynomial)

#SINGLE FITTING OF DATA
params, covs  = curve_fit(spec_function, x_data, y_data, p0=[0,0,0,0],sigma=errors,absolute_sigma=True)
a0, a1,a2,a3 = params #need this bit to unpack the tuple 
errs = np.sqrt(np.diag(covs))
a0_err,a1_err,a2_err,a3_err = errs
#calculated fitted data and find the residuals etc
fit_data = spec_function(x_data, *params)
residuals = y_data - fit_data 
chi_squared = np.sum((residuals / errors) ** 2)
dof = len(y_data) - len(params)  # Degrees of freedom
reduced_chi_squared = chi_squared / dof

print(f"Estimated Single Fit Parameters: \n a0 = {a0}+/-{a0_err}, a1 = {a1}+/-{a1_err}, a2 = {a2}+/-{a2_err}, a3 = {a3}+/-{a3_err} \n rChi2 = {reduced_chi_squared}")

#MONTE CARLO METHOD
N = 10000
a_samples = []
a1_samples = []
a2_samples = []
a3_samples = []
a_errs_mc = []
mc_solutions = [(spec_function(interpolation_range, *params))]

for i in range(N):
    y_mc = y_data +np.random.normal(size=len(y_data), scale=errors )
    try:
        params_mc, covs_mc = curve_fit(spec_function, x_data, y_mc, p0=[a0,a1,a2,a3],sigma=errors, absolute_sigma=True)
        a0_mc, a1_mc,a2_mc,a3_mc = params_mc #need this bit to unpack the tuple 
        errs_mc = np.sqrt(np.diag(covs_mc))
        a0_err_mc,a1_err_mc,a2_err_mc,a3_err_mc = errs_mc
        #calculated fitted data and find the residuals etc
        fit_data_mc = spec_function(interpolation_range, *params_mc)
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
print(f'fractional uncertainty along interpolation range at specified energy is {np.mean(mc_fractional_uncert[interpolation_range_start:interpolation_range_end])}')

#calculated fitted data for montecarlo and find the chi squared for MC
fit_data_mc = spec_function(x_data,  a_mc,a1_mc,a2_mc,a3_mc)
residuals_mc = y_data - fit_data_mc 
chi_squared_mc = np.sum((residuals_mc / errors) ** 2)
reduced_chi_squared_mc = chi_squared_mc / dof
print(f"Estimated MC Parameters: \n a0 = {a_mc}, a1 = {a1_mc}, a2 = {a2_mc}, a3 = {a3_mc} \n rChi2 = {reduced_chi_squared_mc} ")

#plot MC data
plt.scatter(x_data, y_data, label="Data", c='r',marker='o',lw=2)
plt.plot(interpolation_range, mc_solutions_mean, label="Fitted Curve", color='blue')
plt.fill_between(interpolation_range, mc_solutions_mean-mc_solutions_std_dev, mc_solutions_mean+mc_solutions_std_dev, step='post', alpha=0.25)
plt.errorbar(x_data, y_data, yerr=errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
plt.legend()
plt.xlim(0,2000)
plt.ylim(-0.05,0.3)
plt.xlabel("Gamma energy (keV)")
plt.ylabel("Efficiency")
#plt.savefig("mc_function.png")
plt.close()

#plot MC residuals 
plt.scatter(x_data, residuals_mc,c='r',marker='o',lw=2)
plt.errorbar(x_data, residuals_mc,yerr=errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
plt.title("Plot of the residual of the fit")
plt.xlabel("Time (s)")
plt.ylabel("Residual")
#plt.savefig("mc_residuals.png")
plt.close()

# plot distribution of a0 params
plt.hist(a_samples, 100)
plt.title("Plot of parameter distribution")
plt.xlabel("Parameter value")
plt.ylabel("Frequency")
#plt.savefig("mc_param_dist.png")
plt.close()