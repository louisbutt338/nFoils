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

#define efficiency polynomial
def spec_function(energy,a0,a1,a2,a3):
    polynomial = a0 + a1*np.log(energy)**1 + a2*np.log(energy)**2 + a3*np.log(energy)**3 
    return np.exp(polynomial)

#get data
folderpath = "/Users/ljb841@student.bham.ac.uk/gamma_spec/deuteron_hpge/hpge_results_b03_291124/calibration"
filename = "1cm.txt"
def parse_txt(file):
    file_path = f"{folderpath}/{file}"
    with open(file_path,'r') as data_file:
        file_contents = data_file.readlines()       
    #return x,y,errors

#experimental data
x_data = [121.7817
         ,244.6974
         ,344.2785
         ,443.9606
         ,778.9045
         ,867.3800
         ,964.0570
         ,1085.837
         ,1112.076
         ,1408.013]
y_data = [0.135849839 
         ,0.103919441 
         ,0.078887158
         ,0.066050731
         ,0.040737369
         ,0.037480357
         ,0.035006829
         ,0.033901763
         ,0.03159819
         ,0.026876227]
errors = [5.95E-03 
         ,3.80E-03 
         ,3.14E-03
         ,2.55E-03
         ,1.47E-03
         ,1.20E-03
         ,1.28E-03
         ,1.33E-03
         ,1.16E-03
         ,9.43E-04]

# #model data 
# x_data = [100
#          ,150
#          ,200
#          ,400
#          ,600
#          ,800
#          ,1000
#          ,1500]
# y_data = [1.98E-01
#          ,1.63E-01
#          ,1.31E-01
#          ,6.75E-02
#          ,4.64E-02
#          ,3.63E-02
#          ,2.99E-02
#          ,2.14E-02]
# errors = [1.98E-02
#          ,1.63E-02
#          ,1.31E-02
#          ,6.75E-03
#          ,4.64E-03
#          ,3.63E-03
#          ,2.99E-03
#          ,2.14E-03]

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
N = 10
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