import os
import numpy as np # type: ignore
import sys
import statistics
import matplotlib.font_manager # type: ignore
import matplotlib # type: ignore
import matplotlib.pyplot as plt # type: ignore
import json
#from matplotlib.pyplot import cm
#plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Helvetica"
plt.rcParams["font.size"] = 22
plt.rcParams["font.weight"] = "normal"
from collections import Counter
from itertools import count
from random import randint
from pathlib import Path
import math

################# DIRECTORY INPUTS #################

working_directory = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/analysis'
exp_results_directory = '/Users/ljb841@student.bham.ac.uk/gamma_spec/proton_hpge/experimental_activities'

################# FLUX NORMALISATION CALCS ################# 

flux_norm_mean = 0.02756
flux_norm_error = 0.004
calculated_flux_uncertainty_frac = flux_norm_error/flux_norm_mean
n_per_s_norm_mean = 0.0751237
n_per_s_error = 0.0751237* np.sqrt((1147.24/199358.77)**2+(0.004/0.02756)**2+(59186.1/73137.1)**2)

def estimated_10ua_flux(flux):
    return flux*6.24151e13*1.37267e-05
def estimated_10ua_li_n_per_s(flux):
    return flux*6.24151e13*3.28565E-03
print(f"max flux at 10uA proton current: {estimated_10ua_flux(flux_norm_mean):.3e} +- {estimated_10ua_flux(flux_norm_error):.3e} n/cm2/s")
print(f"max n/s from lithium at 10uA proton current: {estimated_10ua_li_n_per_s(n_per_s_norm_mean):.3e} +- {estimated_10ua_li_n_per_s(n_per_s_error):.3e} n/s")

################# LIST OF ANALYSED ISOTOPES ################# 

experiment = 'proton_march24'
model_dataset = 'seperated_foils'

experimental_analysis_lib = 'root'



#####################################

model_results_path = f"{working_directory}/calculated_activities/{experiment}/{model_dataset}.json"
model_results_data = json.load(open(model_results_path))
isotope_list = model_results_data["isotope_data"]["isotope_list"]
foil_weight_normalisation = model_results_data["isotope_data"]["associated_foil_weights"]

calculated_endfb8_activities = model_results_data["endfb8"]["activity"]
calculated_endfb8_uncertainties = model_results_data["endfb8"]["uncertainty"]
calculated_irdff2_activities = model_results_data["irdff2"]["activity"]
calculated_irdff2_uncertainties = model_results_data["irdff2"]["uncertainty"]
calculated_tendl21_activities = model_results_data["tendl21"]["activity"]
calculated_tendl21_uncertainties = model_results_data["tendl21"]["uncertainty"]


exp_results_path =  f"{exp_results_directory}/{experimental_analysis_lib}_activities.json" 
exp_results_data = json.load(open(exp_results_path))
experimental_activities = exp_results_data["activity"]
experimental_uncertainties = exp_results_data["uncertainty"]

################# CALCULATED UNCERTAINTIES FUNCTION ################# 

def total_calculated_uncerts(calc_fispact_uncerts,calc_activities):
    calculated_uncertainties_frac = []
    for i in np.linspace(0,len(calc_activities)-1,len(calc_activities)):
        frac_fispact_uncerts = calc_fispact_uncerts[int(i)]/calc_activities[int(i)]
        squares_calc_uncerts = frac_fispact_uncerts**2 + calculated_flux_uncertainty_frac**2 
        calculated_uncertainties_frac.append(np.sqrt((squares_calc_uncerts)))
    return calculated_uncertainties_frac


################# C/E FUNCTIONS ################# 

def c_over_e(calculated_activities):
    c_over_e_array = []
    for i in np.linspace(0,len(experimental_activities)-1,len(experimental_activities)):
        c_over_e =  (foil_weight_normalisation[int(i)]*calculated_activities[int(i)])/experimental_activities[int(i)]
        c_over_e_array.append(c_over_e)
    return c_over_e_array

def c_over_e_uncerts(calculated_uncertainties_frac):
    c_over_e_uncerts = []
    for i in np.linspace(0,len(experimental_activities)-1,len(experimental_activities)):
        c_error = calculated_uncertainties_frac[int(i)]
        e_error = experimental_uncertainties[int(i)]/experimental_activities[int(i)]

        ce_error =  np.sqrt( c_error**2 + e_error**2 )
        #ce_error =  np.sqrt( c_error**2 + e_error**2 ) * c_over_e[int(i)]
        #ce_error =  ( c_error + e_error ) * c_over_e[int(i)]
    
        c_over_e_uncerts.append(ce_error)
    return c_over_e_uncerts

################# PLOTTING ################# 

new_order = [9,2,11,12, 7,5,6,  8,4,0,3,1,13,10,14]

new_isotope_list = [isotope_list[i] for i in new_order]
new_ce_results    =  [1.03*c_over_e(calculated_tendl21_activities)[i] for i in new_order]
ce_results_irdff  =  [1.03*c_over_e(calculated_irdff2_activities)[i] for i in new_order]
ce_results_endfb8 =  [1.03*c_over_e(calculated_endfb8_activities)[i] for i in new_order]
print(ce_results_endfb8[12])

new_ce_errors =     [c_over_e_uncerts(total_calculated_uncerts(calculated_tendl21_uncertainties,calculated_tendl21_activities))[i] for i in new_order]
ce_errors_irdff =   [c_over_e_uncerts(total_calculated_uncerts(calculated_irdff2_uncertainties, calculated_irdff2_activities ))[i] for i in new_order]
ce_errors_endfb8 =  [c_over_e_uncerts(total_calculated_uncerts(calculated_endfb8_uncertainties, calculated_endfb8_activities ))[i] for i in new_order]

fig, ax1 = plt.subplots()
ax1.set_xlabel('Neutron-transmuted isotopes') 
ax1.set_ylabel('C/E')
ax1.tick_params(axis='y')
ax1.set_xticks(np.arange(0, len(new_order), step=1))

#ax1.set_xscale("log")
ax1.set_ylim(0.1,2)
#ax1.set_ylim(0,4)
ax1.set_yticks([0.2,0.6,1,1.4,1.8])
#ax1.set_yscale("log")

ax1.scatter(new_isotope_list, new_ce_results, s=40 , c='b', linewidth=2,label='TENDL-2021')
ax1.errorbar(new_isotope_list,new_ce_results,new_ce_errors,fmt='none',lw=2,capsize=2,color='Black',zorder=-1)
ax1.set_xlim(-0.5,len(new_order)-0.5)

ax2 = ax1.twiny()
ax2.scatter(new_isotope_list, ce_results_irdff, s=40 , c='magenta', linewidth=2,label='IRDFF-II')
ax2.errorbar(new_isotope_list,ce_results_irdff,ce_errors_irdff,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax2.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax2.set_xlim(-0.3,len(new_order)-0.3)

ax3 = ax1.twiny()
ax3.scatter(new_isotope_list, ce_results_endfb8, s=40 , c='green', linewidth=2,label='ENDF/B-VIII')
ax3.errorbar(new_isotope_list,ce_results_endfb8,ce_errors_endfb8,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax3.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax3.set_xlim(-0.7,len(new_order)-0.7)

ax1.plot([-1,17], np.ones(2), 'Black', ls='--',linewidth=1.5)
ax1.fill_between([-1,17], 1-0.14513788, 1+0.14513788,facecolor='lightgrey',alpha=0.5)

ax1.set_xticklabels(new_isotope_list,rotation=45)
ax1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.90),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
ax2.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
ax3.legend(loc="upper left", bbox_to_anchor=(0.02, 0.82),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
fig.set_size_inches((17, 6))
fig.savefig(os.path.join(f"{working_directory}/ce_plots", 'CE_plot_root_sep_foils_test_int.png'), transparent=False, bbox_inches='tight')

################# WEIGHTED AVES ################# 

def weighted_ce(ce_value_array,ce_error_array):
    weights = []
    weighted_values = []
    for i in np.linspace(0,len(ce_value_array)-1,len(ce_value_array)):
        iterator = int(i)
        weight = 1/((ce_error_array[iterator])**2)
        weighted_value = weight*ce_value_array[iterator]
        weights.append(weight)
        weighted_values.append(weighted_value)
    summed_weights = np.sum(weights)
    summed_weighted_values = np.sum(weighted_values)
    weighted_ce_result = summed_weighted_values/summed_weights
    weighted_ce_error = 1/np.sqrt(summed_weights)
    return weighted_ce_result,weighted_ce_error

print(f"weighted C/E for {new_isotope_list[6:14]} is {weighted_ce(new_ce_results[6:14],new_ce_errors[6:14])[0]} +- {weighted_ce(new_ce_results[6:14],new_ce_errors[6:14])[1]}")

