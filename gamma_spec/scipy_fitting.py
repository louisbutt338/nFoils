#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 12:00:16 2025

@author: ethansumner
"""
import numpy as np # type: ignore
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import pandas as pd 
import scipy.optimize as opt

###########################
####### USER INPUTS #######

#experimental data
b03_endcap_exp_data = {
  121.7817 : {"efficiency" : 0.135849839, "uncertainty": 7.74E-03 },
  244.6974 : {"efficiency" : 0.103919441, "uncertainty": 5.94E-03 },
  344.2785 : {"efficiency" : 0.078887158, "uncertainty": 4.50E-03 },
  443.9606 : {"efficiency" : 0.066050731, "uncertainty": 3.79E-03 },
  778.9045 : {"efficiency" : 0.040737369, "uncertainty": 2.33E-03 },
  867.3800 : {"efficiency" : 0.037480357, "uncertainty": 2.16E-03 },
  964.0570 : {"efficiency" : 0.035006829, "uncertainty": 2.00E-03 },
  1085.837 : {"efficiency" : 0.033901763, "uncertainty": 1.94E-03 },
  1112.076 : {"efficiency" : 0.03159819 , "uncertainty": 1.80E-03 },
  1408.013 : {"efficiency" : 0.026876227, "uncertainty": 1.53E-03 }}
g11_endcap_exp_data = {
  121.7817 : {"efficiency" : 0.165991862, "uncertainty": 9.46E-03 },
  244.6974 : {"efficiency" : 0.089447323, "uncertainty": 5.11E-03 },
  344.2785 : {"efficiency" : 0.054542061, "uncertainty": 3.11E-03 },
  443.9606 : {"efficiency" : 0.041743764, "uncertainty": 2.40E-03 },
  778.9045 : {"efficiency" : 0.019949091, "uncertainty": 1.14E-03 },
  867.3800 : {"efficiency" : 0.017064127, "uncertainty": 9.88E-04 },
  964.0570 : {"efficiency" : 0.015357942, "uncertainty": 8.77E-04 },
  1085.837 : {"efficiency" : 0.015114318, "uncertainty": 8.66E-04 },
  1112.076 : {"efficiency" : 0.012635361, "uncertainty": 7.22E-04 },
  1408.013 : {"efficiency" : 0.010719542, "uncertainty": 6.12E-04 }
}
# model data 
b03_endcap_model_data = {
  100:  {"efficiency" : 1.98E-01, "uncertainty": 1.98E-03 },
  150:  {"efficiency" : 1.63E-01, "uncertainty": 1.63E-03 },
  200:  {"efficiency" : 1.31E-01, "uncertainty": 1.31E-03 },
  400:  {"efficiency" : 6.75E-02, "uncertainty": 6.75E-04 },
  600:  {"efficiency" : 4.64E-02, "uncertainty": 4.64E-04 },
  800:  {"efficiency" : 3.63E-02, "uncertainty": 3.63E-04 },
  1000: {"efficiency" : 2.99E-02, "uncertainty": 2.99E-04 },
  1500: {"efficiency" : 2.14E-02, "uncertainty": 2.14E-04 }}
g11_endcap_model_data = {
  100:  {"efficiency" : 0.20975 , "uncertainty": 0.20975E-01  },
  150:  {"efficiency" : 0.16322 , "uncertainty": 0.16322E-01  },
  200:  {"efficiency" : 0.117447, "uncertainty": 0.117447E-01 },
  400:  {"efficiency" : 0.045982, "uncertainty": 0.045982E-01 },
  600:  {"efficiency" : 0.027943, "uncertainty": 0.027943E-01 },
  800:  {"efficiency" : 0.020153, "uncertainty": 0.020153E-01 },
  1000: {"efficiency" : 0.015819, "uncertainty": 0.015819E-01 },
  1500: {"efficiency" : 0.010543, "uncertainty": 0.010543E-01 }}

# select which dataset
x_data = [i for i in g11_endcap_exp_data.keys()]
y_data = [g11_endcap_exp_data[i]["efficiency" ] for i in x_data]
errors = [g11_endcap_exp_data[i]["uncertainty"] for i in x_data]

###########################
###########################

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

# #plot data 
# plt.scatter(x_data, y_data, label="Data", c='r',marker='o',lw=2)
# plt.plot(np.arange(1,2000), spec_function(np.arange(1,2000), a0,a1,a2,a3), label="Fitted Curve", color='blue')
# plt.errorbar(x_data, y_data, yerr=errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
# plt.legend()
# plt.xlim(1,2000)
# plt.ylim(1e-2,2e-1)
# plt.xlabel("Gamma energy (keV)")
# plt.ylabel("Efficiency")
# plt.savefig("function.png")
# plt.close()
# #plot residuals 
# plt.scatter(x_data, residuals,c='r',marker='o',lw=2)
# plt.errorbar(x_data, residuals,yerr=errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
# plt.title("Plot of the residual of the fit")
# plt.xlabel("Time (s)")
# plt.ylabel("Residual")
# plt.savefig("residuals.png")
# plt.close()

#MONTE CARLO METHOD
N = 100
a_samples = []
a1_samples = []
a2_samples = []
a3_samples = []
a_covs_samples = []


for i in range(N):
    y_mc = y_data +np.random.normal(size=len(y_data), scale=errors )
    try:
        params_mc, covs_mc = curve_fit(spec_function, x_data, y_mc, p0=[a0,a1,a2,a3],sigma=errors, absolute_sigma=True)
        #a_samples, a1_samples,a2_samples,a3_samples = params_mc #need this bit to unpack the tuple 
        a_samples.append(params_mc[0])
        a1_samples.append(params_mc[1])
        a2_samples.append(params_mc[2])
        a3_samples.append(params_mc[3])
        a_covs_samples.append(covs_mc[0][0])

    except RuntimeError:
        pass  # Skip failed fits

# Compute mean parameter values and uncertainties
a_mc =  np.mean(a_samples)
a1_mc = np.mean(a1_samples)
a2_mc = np.mean(a2_samples)
a3_mc = np.mean(a3_samples)
a_err_mc  = np.std(a_samples)/np.sqrt(N)
a1_err_mc = np.std(a1_samples)/np.sqrt(N)
a2_err_mc = np.std(a2_samples)/np.sqrt(N)
a3_err_mc = np.std(a3_samples)/np.sqrt(N)

#calculated fitted data and find the residuals etc for MC
fit_data_mc = spec_function(x_data,  a_mc,a1_mc,a2_mc,a3_mc)
residuals_mc = y_data - fit_data_mc 
chi_squared_mc = np.sum((residuals_mc / errors) ** 2)
reduced_chi_squared_mc = chi_squared_mc / dof

print(f"Estimated MC Parameters: \n a0 = {a_mc}+/-{a_err_mc}, a1 = {a1_mc}+/-{a1_err_mc}, a2 = {a2_mc}+/-{a2_err_mc}, a3 = {a3_mc}+/-{a3_err_mc} \n rChi2 = {reduced_chi_squared_mc} ")

#plot MC data
plt.scatter(x_data, y_data, label="Data", c='r',marker='o',lw=2)
plt.plot(np.arange(1,2000), spec_function(np.arange(1,2000), a_mc,a1_mc,a2_mc,a3_mc), label="Fitted Curve", color='blue')
plt.errorbar(x_data, y_data, yerr=errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
plt.legend()
plt.xlim(0,2000)
#plt.ylim(1e-2,2e-1)
plt.xlabel("Gamma energy (keV)")
plt.ylabel("Efficiency")
plt.savefig("mc_function.png")
plt.close()
#plot MC residuals 
plt.scatter(x_data, residuals_mc,c='r',marker='o',lw=2)
plt.errorbar(x_data, residuals_mc,yerr=errors,lw=2,capsize=2,color='k',zorder=-1,fmt='none')
plt.title("Plot of the residual of the fit")
plt.xlabel("Time (s)")
plt.ylabel("Residual")
plt.savefig("mc_residuals.png")
plt.close()
# plot distribution of params
plt.hist(a_samples, 100)
plt.title("Plot of parameter distribution")
plt.xlabel("Parameter value")
plt.ylabel("Frequency")
plt.savefig("mc_param_dist.png")
plt.close()