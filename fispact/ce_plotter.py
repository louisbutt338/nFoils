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
calculated_results_file = 'calculated_activities_unfolded'
experimental_results_file = 'experimental_activities_apr25'

# select analysis method out of root/interspec
experimental_analysis_method = 'root'

# FLUX NORMALISATION for C results - calculated using be7 calcs from the xs_calculator program
#flux_norm_mean = 0.702246
#flux_norm_error = 0.1496197
#flux_norm_error = 0.5 * 0.040255/0.169795

flux_norm_mean = 1
#flux_norm_mean = 0.6876
flux_percentage_error = 0.25358 #0.25358 ave for unfolded spectrum #0.2197 is flux estimate uncertainty

#####################################
#####################################
#####################################

# extract data for calculated activites
model_results_path = f"{working_directory}/{calculated_results_file}.json"
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
exp_results_path =  f"{working_directory}/{experimental_results_file}.json" 
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
        if experiment == 'proton_march24':
            fispact_error = calc_uncerts[i]/calc_activities[i]
        if experiment == 'deuteron_nov24':
            fispact_error = calc_uncerts[i]
        c_error = np.sqrt( fispact_error**2 + flux_percentage_error**2)
        e_error = experimental_uncertainties[i]/experimental_activities[i]
        ce_error =  np.sqrt( c_error**2 + e_error**2)
        #ce_error =  np.sqrt( c_error**2 + e_error**2 ) * c_over_e[i]
        #ce_error =  ( c_error + e_error ) * c_over_e[i]
        c_over_e_uncerts.append(ce_error)
    return c_over_e_uncerts

#reorder into capture-to-threshold and perform calculations for the foils (top: protons, bottom:deuterons)
new_order = [10,2,12,13, 7,5,6,9 ,4,0,3,1,14,11]
#new_order = [10,15,19,2,8, 16,6,7, 9,12,13,4,5,0, 3, 1,14,17,11,18]
new_isotope_list = [isotope_list_mathmode[i] for i in new_order]

ce_results_tendl  = [(flux_norm_mean)*c_over_e(calculated_tendl21_activities)[i] for i in new_order]
ce_results_irdff  = [(flux_norm_mean)*c_over_e(calculated_irdff2_activities)[i]  for i in new_order]
ce_results_endfb8 = [(flux_norm_mean)*c_over_e(calculated_endfb8_activities)[i]  for i in new_order]
ce_errors_tendl =   [(flux_norm_mean)*c_over_e(calculated_tendl21_activities)[i] *c_over_e_uncerts(calculated_tendl21_uncertainties,calculated_tendl21_activities)[i] for i in new_order]
ce_errors_irdff =   [(flux_norm_mean)*c_over_e(calculated_irdff2_activities)[i]  *c_over_e_uncerts(calculated_irdff2_uncertainties, calculated_irdff2_activities )[i] for i in new_order]
ce_errors_endfb8 =  [(flux_norm_mean)*c_over_e(calculated_endfb8_activities)[i]  *c_over_e_uncerts(calculated_endfb8_uncertainties, calculated_endfb8_activities )[i] for i in new_order]

print(f"{new_isotope_list[0]} C/E is {ce_results_irdff[0]} pm {ce_errors_irdff[0]}")
#for i in new_order:
#    print(f"{new_isotope_list[i]} C/E is {ce_results_endfb8[i]} pm {ce_errors_endfb8[i]}")

# plot c/e diagram with be7 and uncertainty as the error bar
#fig, ax1 = plt.subplots()
plot_split_integer = 6
fig, (ax1,ax4) = plt.subplots(1,2,figsize=(12,6),gridspec_kw={'width_ratios': [len(new_order[:plot_split_integer]), len(new_order[plot_split_integer:])]})
fig.supxlabel('Neutron-transmuted isotopes',x=0.5,y=-0.14) 
fig.supylabel('C/E',x=0.06,y=0.5)
plt.subplots_adjust(wspace=0.05)

# plotting the first half of plot (significant capture reactions)
ax1.tick_params(axis='y',bottom=False,left=True,labelleft=True)
ax1.set_xticks(np.arange(0, len(new_order[:plot_split_integer]), step=1))
ax1.set_ylim(0,25)
#ax1.set_yticks([0,0.5,1,1.5,2,2.5,3])
ax1.scatter (new_isotope_list[:plot_split_integer],ce_results_tendl[:plot_split_integer], s=40 , c='b', linewidth=2,label='TENDL-2021')
ax1.errorbar(new_isotope_list[:plot_split_integer],ce_results_tendl[:plot_split_integer],ce_errors_tendl[:plot_split_integer],fmt='none',lw=2,capsize=2,color='Black',zorder=-1)
ax1.set_xlim(-0.5,len(new_order[:plot_split_integer])-0.5)
ax1.set_xticklabels(new_isotope_list[:plot_split_integer],rotation=45)
ax1.plot([-1,len(new_order)], np.ones(2), 'Black', ls='--',linewidth=1.5)
ax1.fill_between([-1,len(new_order)], 1-flux_percentage_error, 1+flux_percentage_error,facecolor='lightcoral',alpha=0.3)
ax2 = ax1.twiny()
ax2.scatter (new_isotope_list[:plot_split_integer],ce_results_irdff[:plot_split_integer], s=40 , c='magenta', linewidth=2,label='IRDFF-II')
ax2.errorbar(new_isotope_list[:plot_split_integer],ce_results_irdff[:plot_split_integer],ce_errors_irdff[:plot_split_integer],fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax2.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax2.set_xlim(-0.3,len(new_order[:plot_split_integer])-0.3)
ax3 = ax1.twiny()
ax3.scatter (new_isotope_list[:plot_split_integer],ce_results_endfb8[:plot_split_integer], s=40 , c='green', linewidth=2,label='ENDF/B-VIII')
ax3.errorbar(new_isotope_list[:plot_split_integer],ce_results_endfb8[:plot_split_integer],ce_errors_endfb8[:plot_split_integer],fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax3.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax3.set_xlim(-0.7,len(new_order[:plot_split_integer])-0.7)

# plotting the second half of plot (significant threshold reactions)
ax4.tick_params(axis='y',right=True,labelright=True,left=False,labelleft=False,bottom=False)
ax4.set_xticks(np.arange(0, len(new_order[plot_split_integer:]), step=1))
ax4.set_ylim(0.5,3)
#ax1.set_yticks([0,0.5,1,1.5,2,2.5,3])
ax4.scatter (new_isotope_list[plot_split_integer:],ce_results_tendl[plot_split_integer:], s=40 , c='b', linewidth=2,label='TENDL-2021')
ax4.errorbar(new_isotope_list[plot_split_integer:],ce_results_tendl[plot_split_integer:],ce_errors_tendl[plot_split_integer:],fmt='none',lw=2,capsize=2,color='Black',zorder=-1)
ax4.set_xlim(-0.5,len(new_order[plot_split_integer:])-0.5)
ax4.set_xticklabels(new_isotope_list[plot_split_integer:],rotation=45)
ax4.plot([-1,len(new_order[plot_split_integer:])], np.ones(2), 'Black', ls='--',linewidth=1.5)
ax4.fill_between([-1,len(new_order[plot_split_integer:])], 1-flux_percentage_error, 1+flux_percentage_error,facecolor='lightcoral',alpha=0.3)
ax5 = ax4.twiny()
ax5.scatter (new_isotope_list[plot_split_integer:],ce_results_irdff[plot_split_integer:], s=40 , c='magenta', linewidth=2,label='IRDFF-II')
ax5.errorbar(new_isotope_list[plot_split_integer:],ce_results_irdff[plot_split_integer:],ce_errors_irdff[plot_split_integer:],fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax5.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax5.set_xlim(-0.3,len(new_order[plot_split_integer:])-0.3)
ax6 = ax4.twiny()
ax6.scatter (new_isotope_list[plot_split_integer:],ce_results_endfb8[plot_split_integer:], s=40 , c='green', linewidth=2,label='ENDF/B-VIII')
ax6.errorbar(new_isotope_list[plot_split_integer:],ce_results_endfb8[plot_split_integer:],ce_errors_endfb8[plot_split_integer:],fmt='none',lw=2,capsize=2,color='black',zorder=-1)
ax6.tick_params(top=False, labeltop=False, bottom=True, labelbottom=False)
ax6.set_xlim(-0.7,len(new_order[plot_split_integer:])-0.7)

# legend and saving figure
ax1.legend(loc="upper left", bbox_to_anchor=(0.62, 0.90),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
ax2.legend(loc="upper left", bbox_to_anchor=(0.62, 0.98),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
ax3.legend(loc="upper left", bbox_to_anchor=(0.62, 0.82),handlelength=0,borderaxespad=0, frameon=False,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
fig.set_size_inches((17, 6))
fig.savefig(os.path.join(f"{working_directory}/ce_plots", f'{experimental_analysis_method}_unfolded_apr2025_test.png'), transparent=False, bbox_inches='tight')

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
print(f"weighted C/E for {new_isotope_list[:]} is {weighted_ce(ce_results_tendl[:],ce_errors_tendl[:])[0]} +- {weighted_ce(ce_results_tendl[:],ce_errors_tendl[:])[1]}")

