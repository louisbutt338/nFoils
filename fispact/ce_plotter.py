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

#####################################
################# INPUTS ############
#####################################

# select working directory with 'calculated_activities' and 'experimental_activities' inside
experiment = 'proton_march24'
working_directory = f'/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/analysis/{experiment}'

# select analysis method out of root/interspec
experimental_analysis_method = 'root'


#####################################
#####################################
#####################################

# extract data for calculated activites
model_results_path = f"{working_directory}/calculated_activities.json"
model_results_data = json.load(open(model_results_path))
isotope_list = model_results_data.keys()
foil_weight_normalisation = [model_results_data[key]["foil_weight"] for key in isotope_list]
isotope_list_mathmode = [model_results_data[key]["mathmode_name"] for key in isotope_list]
calculated_endfb8_activities = [model_results_data[key]["endfb8_values"][0] for key in isotope_list]
calculated_endfb8_uncertainties = [model_results_data[key]["endfb8_values"][1] for key in isotope_list]
calculated_irdff2_activities = [model_results_data[key]["irdff2_values"][0] for key in isotope_list]
calculated_irdff2_uncertainties = [model_results_data[key]["irdff2_values"][1] for key in isotope_list]
calculated_tendl21_activities = [model_results_data[key]["tendl21_values"][0] for key in isotope_list]
calculated_tendl21_uncertainties = [model_results_data[key]["tendl21_values"][1] for key in isotope_list]

# extract data for experimental activities
exp_results_path =  f"{working_directory}/experimental_activities.json" 
exp_results_data = json.load(open(exp_results_path))
experimental_activities = [exp_results_data[key][f"{experimental_analysis_method}_values"][0] for key in isotope_list]
experimental_uncertainties = [exp_results_data[key][f"{experimental_analysis_method}_values"][1] for key in isotope_list]

#c/e function
def c_over_e(calculated_activities):
    c_over_e_array = []
    for i in np.arange(len(isotope_list)):
        c_over_e =  (foil_weight_normalisation[i]*calculated_activities[i])/experimental_activities[i]
        c_over_e_array.append(c_over_e)
    return c_over_e_array

#c/e uncertainity function
def c_over_e_uncerts(calc_uncerts,calc_activities):
    c_over_e_uncerts = []
    for i in np.arange(len(isotope_list)):
        #c_error = calculated_uncertainties_frac[i]
        c_error = calc_uncerts[i]/calc_activities[i]
        e_error = experimental_uncertainties[i]/experimental_activities[i]
        ce_error =  np.sqrt( c_error**2 + e_error**2 )
        #ce_error =  np.sqrt( c_error**2 + e_error**2 ) * c_over_e[i]
        #ce_error =  ( c_error + e_error ) * c_over_e[i]
        c_over_e_uncerts.append(ce_error)
    return c_over_e_uncerts

# FLUX NORMALISATION CALCS with be7 result - leave n/s calculations to the be7 xs_calculator program

be7_norm_factor = model_results_data["Be7"]["endfb8_values"][0] / exp_results_data["Be7"][f"{experimental_analysis_method}_values"][0]
be7_frac_uncertainty = exp_results_data["Be7"][f"{experimental_analysis_method}_values"][1] / exp_results_data["Be7"][f"{experimental_analysis_method}_values"][0]
flux_norm_mean = be7_norm_factor
flux_norm_error = (be7_frac_uncertainty*be7_norm_factor)
def estimated_10ua_flux(flux_norm):
    if experiment == 'proton_march24':
        flux_10ua = 0.02756*flux_norm*6.24151e13*1.37267e-05
    if experiment == 'deuteron_nov24':
        flux_10ua = flux_norm*3.42800e9
    return flux_10ua
print(f"max flux on Fe foil at 10uA proton current: {estimated_10ua_flux(flux_norm_mean):.3e} +- {estimated_10ua_flux(flux_norm_error):.3e} n/cm2/s")

#reorder into capture-to-threshold and perform calculations for the foils
new_order = [10,2,12,13, 7,5,6,9 ,4,0,3,1,14,11]
new_isotope_list = [isotope_list_mathmode[i] for i in new_order]
ce_results_tendl  = [(1/be7_norm_factor)*c_over_e(calculated_tendl21_activities)[i] for i in new_order]
ce_results_irdff  = [(1/be7_norm_factor)*c_over_e(calculated_irdff2_activities)[i] for i in new_order]
ce_results_endfb8 = [(1/be7_norm_factor)*c_over_e(calculated_endfb8_activities)[i] for i in new_order]
ce_errors_tendl =   [(1/be7_norm_factor)*c_over_e_uncerts(calculated_tendl21_uncertainties,calculated_tendl21_activities)[i] for i in new_order]
ce_errors_irdff =   [(1/be7_norm_factor)*c_over_e_uncerts(calculated_irdff2_uncertainties, calculated_irdff2_activities )[i] for i in new_order]
ce_errors_endfb8 =  [(1/be7_norm_factor)*c_over_e_uncerts(calculated_endfb8_uncertainties, calculated_endfb8_activities )[i] for i in new_order]


# plot c/e diagram with be7 and uncertainty as the error bar
fig, ax1 = plt.subplots()
ax1.set_xlabel('Neutron-transmuted isotopes') 
ax1.set_ylabel('C/E')
ax1.tick_params(axis='y')
ax1.set_xticks(np.arange(0, len(new_order), step=1))
#ax1.set_xscale("log")
ax1.set_ylim(0.1,2.5)
#ax1.set_ylim(0,4)
ax1.set_yticks([0.2,0.6,1,1.4,1.8,2.2])
#ax1.set_yscale("log")
####
ax1.scatter(new_isotope_list, ce_results_tendl, s=40 , c='b', linewidth=2,label='TENDL-2021')
ax1.errorbar(new_isotope_list,ce_results_tendl,ce_errors_tendl,fmt='none',lw=2,capsize=2,color='Black',zorder=-1)
ax1.set_xlim(-0.5,len(new_order)-0.5)
####
ax2 = ax1.twiny()
ax2.scatter(new_isotope_list, ce_results_irdff, s=40 , c='magenta', linewidth=2,label='IRDFF-II')
ax2.errorbar(new_isotope_list,ce_results_irdff,ce_errors_irdff,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax2.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax2.set_xlim(-0.3,len(new_order)-0.3)
#####
ax3 = ax1.twiny()
ax3.scatter(new_isotope_list, ce_results_endfb8, s=40 , c='green', linewidth=2,label='ENDF/B-VIII')
ax3.errorbar(new_isotope_list,ce_results_endfb8,ce_errors_endfb8,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax3.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax3.set_xlim(-0.7,len(new_order)-0.7)
#####
ax1.plot([-1,17], np.ones(2), 'Black', ls='--',linewidth=1.5)
ax1.fill_between([-1,17], 1-be7_frac_uncertainty, 1+be7_frac_uncertainty,facecolor='lightcoral',alpha=0.3)
#####
ax1.set_xticklabels(new_isotope_list,rotation=45)
ax1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.90),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
ax2.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
ax3.legend(loc="upper left", bbox_to_anchor=(0.02, 0.82),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
fig.set_size_inches((17, 6))
fig.savefig(os.path.join(f"{working_directory}/ce_plots", 'root_sep_foils.png'), transparent=False, bbox_inches='tight')

# calculate weighted averages
def weighted_ce(ce_value_array,ce_error_array):
    weights = []
    weighted_values = []
    for i in np.arange(len(ce_value_array)):
        weight = 1/((ce_error_array[i])**2)
        weighted_value = weight*ce_value_array[i]
        weights.append(weight)
        weighted_values.append(weighted_value)
    summed_weights = np.sum(weights)
    summed_weighted_values = np.sum(weighted_values)
    weighted_ce_result = summed_weighted_values/summed_weights
    weighted_ce_error = 1/np.sqrt(summed_weights)
    return weighted_ce_result,weighted_ce_error
print(f"weighted C/E for {new_isotope_list[6:14]} is {weighted_ce(ce_results_tendl[6:14],ce_errors_tendl[6:14])[0]} +- {weighted_ce(ce_results_tendl[6:14],ce_errors_tendl[6:14])[1]}")

