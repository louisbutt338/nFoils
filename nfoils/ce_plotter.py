""" module for performing ce analysis for some isotopes
"""

import os
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},weight='normal',size=20)
import json
from pathlib import Path
import math

class CEPlotter:
    """ class for doing the ce analysis
    """
    def __init__(self, experiment,calculated_results_file,
                experimental_results_file,plotname,y_upper,
                y_lower,flux_norm_mean,flux_percentage_error,
                first_we,last_we):
        """ Initialise class
        """

        # set params
        self.experiment_directory = experiment
        self.calculated_results_file =  calculated_results_file
        self.experimental_results_file = experimental_results_file
        self.plotname = plotname
        self.flux_norm_mean = flux_norm_mean
        self.flux_percentage_error = flux_percentage_error
        self.first_we = first_we
        self.last_we = last_we
        self.y_axis_upper = y_upper
        self.y_axis_lower = y_lower

    def _c_over_e(self,calc_activities,exp_activities,isotope_list,foil_weight,ssf):
        """ c/e calculation

        Parameters
        ----------
        calc_activities : list
            List of calculated activities
        exp_activities : list
            List of measured activities
        isotope_list : list
            List of the isotopes
        foil_weight : list
            List of the weights of the respective foils
        ssf : list
            List of the self shielding factors for the respective isotopes
        """
        c_over_e_array = []
        for i in np.arange(len(isotope_list)):
            c_over_e =  (foil_weight[i]*ssf[i]*calc_activities[i])/exp_activities[i]
            c_over_e_array.append(c_over_e)
        return c_over_e_array

    def _c_over_e_uncerts(self,calc_uncerts,calc_activities,experimental_activities,
                          experimental_uncertainties,isotope_list,
                          isotopic_spectrum_percentage_uncerts):
        """ c/e uncertainty calculation

        Parameters
        ----------
        calc_uncerts : list
            List of calculated uncertainties
        calc_activities : list
            List of calculated activities
        exp_activities : list
            List of measured activities
        exp_uncertainties : list
            List of measured uncerts
        isotope_list : list
            List of the isotopes
        isotopic_spectrum_percentage_uncerts : list
            List of the uncerts in the spectrum for all the isotopes
        """
        c_over_e_uncerts = []
        for i in np.arange(len(isotope_list)):
            if self.experiment_directory == 'proton_march24':
                fispact_error = calc_uncerts[i]/calc_activities[i]
            if self.experiment_directory == 'deuteron_nov24':
                fispact_error = calc_uncerts[i]
            c_error = np.sqrt( fispact_error**2 
                              + isotopic_spectrum_percentage_uncerts[i]**2)
            e_error = experimental_uncertainties[i]/experimental_activities[i]
            ce_error =  np.sqrt( c_error**2 + e_error**2)
            #ce_error =  np.sqrt( c_error**2 + e_error**2 ) * c_over_e[i]
            #ce_error =  ( c_error + e_error ) * c_over_e[i]
            c_over_e_uncerts.append(ce_error)
        return c_over_e_uncerts

    def _plotter(self,new_order,new_isotope_list,ce_results_1,
                 ce_errors_1,ce_results_2,ce_errors_2,
                 ce_results_3,ce_errors_3,
                 isotopic_flux_percentage_uncerts,library_labels):
        """ plot c/e diagram with flux+spectrum uncertainty as the error bar
        can do three nuclear data libraries

        Parameters
        ----------
        new_order : list
            New order for the isotopes to be plotted in
        new_isotope_list : list
            Isotopes listed in their new order
        ce_results_1 : list
            List of ce results for library 1
        ce_errors_1 : list
            List of ce errors for library 1
        ce_results_2 : list
            List of ce results for library 2
        ce_errors_2 : list
            List of ce errors for library 2
        ce_results_3 : list
            List of ce results for library 3
        ce_errors_3 : list
            List of ce errors for library 3
        isotopic_flux_percentage_uncerts : list
            List of the 
        isotopic_spectrum_percentage_uncerts : list
            List of total spectrum+flux uncerts for each isotope
        library_labels : list
            List of labels of the libraries for the plot
        """
        #initial plotting settings
        plot_split_integer = 5
        fig, (ax1,ax4) = plt.subplots(
            1,2,figsize=(12,6)
            ,gridspec_kw={'width_ratios': [len(new_order[:plot_split_integer])
                                           , len(new_order[plot_split_integer:])]})
        fig.supxlabel('Neutron-transmuted isotopes',x=0.5,y=-0.14) 
        fig.supylabel('C/E',x=0.06,y=0.5)
        plt.subplots_adjust(wspace=0.05)

        # plotting the first half of plot (significant capture reactions)
        ax1.tick_params(axis='y',bottom=False,left=True,labelleft=True,top=True)
        ax1.set_xticks(np.arange(len(new_order[:plot_split_integer]))
                       , labels=new_isotope_list[:plot_split_integer] ,rotation=45)
        ax1.set_ylim(self.y_axis_lower,self.y_axis_upper)
        ax1.scatter (new_isotope_list[:plot_split_integer]
                     ,ce_results_1[:plot_split_integer]
                     , s=40 , c='b', linewidth=2,label=library_labels[0])
        ax1.errorbar(new_isotope_list[:plot_split_integer]
                     ,ce_results_1[:plot_split_integer]
                     ,ce_errors_1[:plot_split_integer]
                     ,fmt='none',lw=2,capsize=2,color='Black',zorder=-1)
        ax1.set_xlim(-0.4,len(new_order[:plot_split_integer])-0.6)
        ax2 = ax1.twiny()
        ax2.scatter (new_isotope_list[:plot_split_integer]
                     ,ce_results_2[:plot_split_integer]
                     , s=40 , c='magenta', linewidth=2,label=library_labels[1])
        ax2.errorbar(new_isotope_list[:plot_split_integer]
                     ,ce_results_2[:plot_split_integer]
                     ,ce_errors_2[:plot_split_integer]
                     ,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
        ax2.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
        ax2.set_xlim(-0.2,len(new_order[:plot_split_integer])-0.4)
        ax3 = ax1.twiny()
        ax3.scatter (new_isotope_list[:plot_split_integer]
                     ,ce_results_3[:plot_split_integer]
                     , s=40 , c='green', linewidth=2,label=library_labels[2])
        ax3.errorbar(new_isotope_list[:plot_split_integer]
                     ,ce_results_3[:plot_split_integer]
                     ,ce_errors_3[:plot_split_integer]
                     ,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
        ax3.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
        ax3.set_xlim(-0.6,len(new_order[:plot_split_integer])-0.8)
        ax1.plot([-1,len(new_order)], np.ones(2), 'Black', ls='--',linewidth=1.5)
        ax1.fill_between(range(len(new_order[:plot_split_integer]))
                         ,[1-x for x in isotopic_flux_percentage_uncerts
                           [:plot_split_integer]]
                         ,[1+x for x in isotopic_flux_percentage_uncerts
                           [:plot_split_integer]]
                         ,facecolor='lightcoral',alpha=0.3,step='mid')

        # plotting the second half of plot (significant threshold reactions)
        ax4.tick_params(axis='y',right=False,labelright=False
                        ,left=False,labelleft=False,bottom=False)
        ax4.set_xticks(np.arange(len(new_order[plot_split_integer:]))
                       ,labels=new_isotope_list[plot_split_integer:],rotation=45)
        ax4.set_ylim(self.y_axis_lower,self.y_axis_upper)
        ax4.scatter (new_isotope_list[plot_split_integer:]
                     ,ce_results_1[plot_split_integer:]
                     , s=40 , c='b', linewidth=2,label=library_labels[0])
        ax4.errorbar(new_isotope_list[plot_split_integer:]
                     ,ce_results_1[plot_split_integer:]
                     ,ce_errors_1[plot_split_integer:]
                     ,fmt='none',lw=2,capsize=2,color='Black',zorder=-1)
        ax4.set_xlim(-0.4,len(new_order[plot_split_integer:])-0.6)
        ax4.set_xticklabels(new_isotope_list[plot_split_integer:],rotation=45)
        ax5 = ax4.twiny()
        ax5.scatter (new_isotope_list[plot_split_integer:]
                     ,ce_results_2[plot_split_integer:]
                     , s=40 , c='magenta', linewidth=2,label=library_labels[1])
        ax5.errorbar(new_isotope_list[plot_split_integer:]
                     ,ce_results_2[plot_split_integer:]
                     ,ce_errors_2[plot_split_integer:]
                     ,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
        ax5.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
        ax5.set_xlim(-0.2,len(new_order[plot_split_integer:])-0.4)
        ax6 = ax4.twiny()
        ax6.scatter (new_isotope_list[plot_split_integer:]
                     ,ce_results_3[plot_split_integer:]
                     , s=40 , c='green', linewidth=2,label=library_labels[2])
        ax6.errorbar(new_isotope_list[plot_split_integer:]
                     ,ce_results_3[plot_split_integer:]
                     ,ce_errors_3[plot_split_integer:]
                     ,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
        ax6.tick_params(top=False, labeltop=False, bottom=False, labelbottom=False)
        ax6.set_xlim(-0.6,len(new_order[plot_split_integer:])-0.8)
        ax4.plot([-1,len(new_order[plot_split_integer:])]
                 , np.ones(2), 'Black', ls='--',linewidth=1.5)
        ax4.fill_between(range(len(new_order[plot_split_integer:]))
                         ,[1-x for x in isotopic_flux_percentage_uncerts
                           [plot_split_integer:]]
                         ,[1+x for x in isotopic_flux_percentage_uncerts
                           [plot_split_integer:]]
                         ,facecolor='lightcoral',alpha=0.3,step='mid')

        # legend and saving figure
        ax4.set_zorder(-1)
        ax1.legend(loc="upper left", bbox_to_anchor=(1.12, 0.90)
                   ,handlelength=0,borderaxespad=0, frameon=False
                   ,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
        ax2.legend(loc="upper left", bbox_to_anchor=(1.12, 0.98)
                   ,handlelength=0,borderaxespad=0, frameon=False
                   ,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
        ax3.legend(loc="upper left", bbox_to_anchor=(1.12, 0.82)
                   ,handlelength=0,borderaxespad=0, frameon=False
                   ,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
        fig.set_size_inches((17, 6))
        fig.savefig(os.path.join(
            f"{self.experiment_directory}", f'{self.plotname}.png')
            ,transparent=False, bbox_inches='tight')

    def _weighted_ce(self,ce_value_array,ce_error_array):
        """ calculate some weighted c/e averages from first_we to last_we

        Parameters
        ----------
        ce_value_array : array
            Array of C/ values to do the maffs on
        ce_error_array : list
            See above but for errors
        """
        weights = []
        weighted_values = []
        for i in np.arange(len(ce_value_array)):
            if ce_value_array[i] <10:
                weight = 1/((ce_error_array[i])**2)
                weighted_value = weight*ce_value_array[i]
                weights.append(weight)
                weighted_values.append(weighted_value)
        summed_weights = np.sum(weights)
        summed_weighted_values = np.sum(weighted_values)
        weighted_ce_result = summed_weighted_values/summed_weights
        weighted_ce_error = 1/np.sqrt(summed_weights)
        return weighted_ce_result,weighted_ce_error
    
    def run(self):
        """ run the analysis
        This is specific to TENDL2021/ENDFB/IRDFF2 so names will need changing 
        if you are using different libraries
        """
        #extract data for calculated activities
        model_results_path = f"{
            self.experiment_directory}/{self.calculated_results_file}.json"
        model_results_data = json.load(open(model_results_path))
        isotope_list = model_results_data.keys()
        foil_weight = [model_results_data[key]["foil_weight"] for key in isotope_list]
        ssf = [model_results_data[key]["ssf_correction"] for key in isotope_list]
        isotope_list_mathmode = [model_results_data[key]["mathmode_name"]
                                  for key in isotope_list]
        calculated_endfb8_a = [model_results_data[key]["endfb8_values"][0]
                                for key in isotope_list]
        calculated_endfb8_u = [model_results_data[key]["endfb8_values"][1]
                                for key in isotope_list]
        calculated_irdff2_a = [model_results_data[key]["irdff2_values"][0]
                                for key in isotope_list]
        calculated_irdff2_u = [model_results_data[key]["irdff2_values"][1]
                                for key in isotope_list]
        calculated_tendl21_a = [model_results_data[key]["tendl21_values"][0]
                                 for key in isotope_list]
        calculated_tendl21_u = [model_results_data[key]["tendl21_values"][1]
                                 for key in isotope_list]

        # calculate spectrum + flux uncertainties for each isotope
        isotopic_spectrum_u = [model_results_data[key]["spectrum_percent_uncert"]
                                for key in isotope_list]
        isotopic_spectrum_u = [np.sqrt(self.flux_percentage_error**2 + i**2)
                                for i in isotopic_spectrum_u]

        #extract data for experimental activities (using average activities)
        exp_results_path =  f"{
            self.experiment_directory}/{self.experimental_results_file}.json" 
        exp_results_data = json.load(open(exp_results_path))
        experimental_a = [np.mean(exp_results_data[key][f"activities"])
                           for key in isotope_list]
        experimental_u = [np.mean(exp_results_data[key][f"activity_uncertainties"])
                           for key in isotope_list]
        
        #reorder results into capture-to-threshold  
        if self.experiment_directory == 'proton_march24':
            new_order = [10,2,12,13, 7,5,6,9 ,4,3,0,1,14,11]
        if self.experiment_directory == 'deuteron_nov24':
            new_order = [10,18,2,15,16, 6,19,7,9,12,13,4,5,0,3,1,17,11]
        new_isotope_list = [isotope_list_mathmode[i] for i in new_order]

        #perform C/E calculations for the foils
        ce_results_tendl  = [
            (self.flux_norm_mean)*
            self._c_over_e(calculated_tendl21_a,experimental_a,isotope_list,foil_weight,ssf)
            [i] for i in new_order]
        ce_results_irdff  = [
            (self.flux_norm_mean)*
            self._c_over_e(calculated_irdff2_a,experimental_a,isotope_list,foil_weight,ssf)
            [i]  for i in new_order]
        ce_results_endfb8 = [
            (self.flux_norm_mean)*
            self._c_over_e(calculated_endfb8_a ,experimental_a,isotope_list,foil_weight,ssf)
            [i]  for i in new_order]
        ce_errors_tendl =   [
            (self.flux_norm_mean)*
            self._c_over_e(calculated_tendl21_a,experimental_a,isotope_list,foil_weight,ssf)[i]
            *self._c_over_e_uncerts(calculated_tendl21_u,calculated_tendl21_a,experimental_a
                                    ,experimental_u,isotope_list,isotopic_spectrum_u)[i]
                                      for i in new_order]
        ce_errors_irdff =   [
            (self.flux_norm_mean)*
            self._c_over_e(calculated_irdff2_a ,experimental_a,isotope_list,foil_weight,ssf)[i]
            *self._c_over_e_uncerts(calculated_irdff2_u, calculated_irdff2_a,experimental_a
                                    ,experimental_u,isotope_list,isotopic_spectrum_u)[i]
                                      for i in new_order]
        ce_errors_endfb8 =  [
            (self.flux_norm_mean)*
            self._c_over_e(calculated_endfb8_a ,experimental_a,isotope_list,foil_weight,ssf)[i]
            *self._c_over_e_uncerts(calculated_endfb8_u, calculated_endfb8_a,experimental_a
                                    ,experimental_u,isotope_list,isotopic_spectrum_u)[i]
                                      for i in new_order]
        isotopic_spectrum_u = [isotopic_spectrum_u[i] for i in new_order]
        
        #print some results if you want
        for i in range(len(new_isotope_list)):
            print(f'*********{new_isotope_list[i]} C/E results')
            print(f"TENDL value is {ce_results_tendl[i]} pm {ce_errors_tendl[i] }")
            print(f"ENDFB value is {ce_results_endfb8[i]} pm {ce_errors_endfb8[i] }")
            print(f"IRDFF value is {ce_results_irdff[i]} pm {ce_errors_irdff[i] }")
        
        # plot results
        library_labels = ['TENDL-2021','IRDFF-II','ENDF/B-VIII']
        self._plotter(new_order,new_isotope_list,ce_results_tendl,
                    ce_errors_tendl,ce_results_irdff,ce_errors_irdff,
                    ce_results_endfb8,ce_errors_endfb8,
                    isotopic_spectrum_u,library_labels)
        
        # do weighted ave calcs and print
        weighted_ave_value = self._weighted_ce(
            ce_results_endfb8[self.first_we:self.last_we],
            ce_errors_endfb8[self.first_we:self.last_we])[0]
        weighted_ave_uncert = self._weighted_ce(
            ce_results_endfb8[self.first_we:self.last_we],
            ce_errors_endfb8[self.first_we:self.last_we])[1]
        print(f"weighted C/E for {new_isotope_list[self.first_we:self.last_we]} is "
              f"{weighted_ave_value} +- {weighted_ave_uncert}")