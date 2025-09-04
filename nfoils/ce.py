""" 
module for performing ce analysis for some isotopes
"""

import json
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=20)


class CEPlotter:
    """ class for doing the ce analysis
    """
    def __init__(self, plotname):
        """ Initialise class
        """

        # set params
        self.plotname = plotname

    def _c_over_e(self,calc_activities,exp_activities,
                  isotope_list,foil_weight,ssf):
        """ c/e calculation from lists of calculated
        and experimental activities

        Parameters
        ----------
        calc_activities : list[float]
            List of calculated activities
        exp_activities : list[float]
            List of measured activities
        isotope_list : list[str]
            List of the isotopes not in mathmode
        foil_weight : list[float]
            List of the weights of the respective foils
        ssf : list[float]
            List of the self shielding factors
            for the respective isotopes

        Returns
        -------
        c_over_e_list : list[float]
            list of c/e values
        """
        c_over_e_list = []
        for i in np.arange(len(isotope_list)):
            c_over_e =  ((foil_weight[i]*ssf[i]*calc_activities[i])/
                         exp_activities[i])
            c_over_e_list.append(c_over_e)
        return c_over_e_list

    def _c_over_e_uncerts(self,calc_uncerts,exp_activities,
                          exp_uncerts,isotope_list,spectrum_flux_frac_uncerts):
        """ c/e uncertainty calculation

        Parameters
        ----------
        calc_uncerts : list[float]
            List of calculated uncertainties
        exp_activities : list[float]
            List of measured activities
        exp_uncerts : list[float]
            List of measured uncerts
        isotope_list : list[str]
            List of the isotopes not in mathmode
        spectrum_flux_frac_uncerts : list[float]
            List of the uncerts in the spectrum for all the isotopes

        Returns
        -------
        c_over_e_uncerts : list[float]
            list of c/e uncertainties
        """
        c_over_e_uncerts = []
        for i in np.arange(len(isotope_list)):
            fispact_error = calc_uncerts[i]
            c_error = np.sqrt(fispact_error**2 
                              +spectrum_flux_frac_uncerts[i]**2)
            e_error = exp_uncerts[i]/exp_activities[i]
            ce_error =  np.sqrt( c_error**2 + e_error**2)
            #ce_error =  np.sqrt( c_error**2 + e_error**2 ) * c_over_e[i]
            #ce_error =  ( c_error + e_error ) * c_over_e[i]
            c_over_e_uncerts.append(ce_error)
        return c_over_e_uncerts

    def _plotter(self,new_order,new_isotope_list,ce_results_1,
                 ce_errors_1,ce_results_2,ce_errors_2,
                 ce_results_3,ce_errors_3,
                 spectrum_flux_frac_uncerts,library_labels,
                 plot_split_integer,y_axis):
        """ plot c/e diagram with flux+spectrum uncertainty as the error bar
        can do three nuclear data libraries

        Parameters
        ----------
        new_order : list[int]
            New order for the isotopes to be plotted in
        new_isotope_list : list[str]
            Isotopes mathmode names listed in their new order
        ce_results_1 : list[float]
            List of ce results for library 1
        ce_errors_1 : list[float]
            List of ce errors for library 1
        ce_results_2 : list[float]
            List of ce results for library 2
        ce_errors_2 : list[float]
            List of ce errors for library 2
        ce_results_3 : list[float]
            List of ce results for library 3
        ce_errors_3 : list[float]
            List of ce errors for library 3
        spectrum_flux_frac_uncerts : list[float]
            List of the total spectrum+flux uncerts for each isotope
        library_labels : list[str]
            List of labels of the libraries for the plot
        plot_split_integer : int
            int to split the plot into thermal/threshold
        y_axis : list[float]
            [lower y value, upper y value] for plot
        """
        #initial plotting settings
        fig, (ax1,ax4) = plt.subplots(
            1,2,figsize=(12,6),
            gridspec_kw={'width_ratios':[len(new_order[:plot_split_integer])
                                         ,len(new_order[plot_split_integer:]
                                              )]})
        fig.supxlabel('Neutron-transmuted isotopes',x=0.5,y=-0.14) 
        fig.supylabel('C/E',x=0.06,y=0.5)
        plt.subplots_adjust(wspace=0.05)

        # plotting the first half of plot (significant capture reactions)
        ax1.tick_params(axis='y',bottom=False,left=True,labelleft=True,
                        top=True)
        ax1.set_xticks(np.arange(len(new_order[:plot_split_integer]))
                       ,labels=new_isotope_list[:plot_split_integer],
                       rotation=45)
        ax1.set_ylim(y_axis[0],y_axis[1])
        ax1.scatter (new_isotope_list[:plot_split_integer]
                     ,ce_results_1[:plot_split_integer]
                     ,s=40 , c='b', linewidth=2,label=library_labels[0])
        ax1.errorbar(new_isotope_list[:plot_split_integer]
                     ,ce_results_1[:plot_split_integer]
                     ,ce_errors_1[:plot_split_integer]
                     ,fmt='none',lw=2,capsize=2,color='Black',zorder=-1)
        ax1.set_xlim(-0.4,len(new_order[:plot_split_integer])-0.6)
        ax2 = ax1.twiny()
        ax2.scatter (new_isotope_list[:plot_split_integer]
                     ,ce_results_2[:plot_split_integer]
                     ,s=40 , c='magenta', linewidth=2,label=library_labels[1])
        ax2.errorbar(new_isotope_list[:plot_split_integer]
                     ,ce_results_2[:plot_split_integer]
                     ,ce_errors_2[:plot_split_integer]
                     ,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
        ax2.tick_params(top=False, labeltop=False, bottom=False,
                        labelbottom=False)
        ax2.set_xlim(-0.2,len(new_order[:plot_split_integer])-0.4)
        ax3 = ax1.twiny()
        ax3.scatter (new_isotope_list[:plot_split_integer]
                     ,ce_results_3[:plot_split_integer]
                     ,s=40 , c='green', linewidth=2,label=library_labels[2])
        ax3.errorbar(new_isotope_list[:plot_split_integer]
                     ,ce_results_3[:plot_split_integer]
                     ,ce_errors_3[:plot_split_integer]
                     ,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
        ax3.tick_params(top=False, labeltop=False, bottom=False,
                        labelbottom=False)
        ax3.set_xlim(-0.6,len(new_order[:plot_split_integer])-0.8)
        ax1.plot([-1,len(new_order)], np.ones(2), 'Black', ls='--',
                 linewidth=1.5)
        ax1.fill_between(range(len(new_order[:plot_split_integer]))
                         ,[1-x for x in spectrum_flux_frac_uncerts
                           [:plot_split_integer]]
                         ,[1+x for x in spectrum_flux_frac_uncerts
                           [:plot_split_integer]]
                         ,facecolor='lightcoral',alpha=0.3,step='mid')

        # plotting the second half of plot (significant threshold reactions)
        ax4.tick_params(axis='y',right=False,labelright=False
                        ,left=False,labelleft=False,bottom=False)
        ax4.set_xticks(np.arange(len(new_order[plot_split_integer:]))
                       ,labels=new_isotope_list[plot_split_integer:],
                       rotation=45)
        ax4.set_ylim(y_axis[0],y_axis[1])
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
        ax5.tick_params(top=False, labeltop=False, bottom=False,
                        labelbottom=False)
        ax5.set_xlim(-0.2,len(new_order[plot_split_integer:])-0.4)
        ax6 = ax4.twiny()
        ax6.scatter (new_isotope_list[plot_split_integer:]
                     ,ce_results_3[plot_split_integer:]
                     , s=40 , c='green', linewidth=2,label=library_labels[2])
        ax6.errorbar(new_isotope_list[plot_split_integer:]
                     ,ce_results_3[plot_split_integer:]
                     ,ce_errors_3[plot_split_integer:]
                     ,fmt='none',lw=2,capsize=2,color='black',zorder=-1)
        ax6.tick_params(top=False, labeltop=False, bottom=False,
                        labelbottom=False)
        ax6.set_xlim(-0.6,len(new_order[plot_split_integer:])-0.8)
        ax4.plot([-1,len(new_order[plot_split_integer:])]
                 , np.ones(2), 'Black', ls='--',linewidth=1.5)
        ax4.fill_between(range(len(new_order[plot_split_integer:]))
                         ,[1-x for x in spectrum_flux_frac_uncerts
                           [plot_split_integer:]]
                         ,[1+x for x in spectrum_flux_frac_uncerts
                           [plot_split_integer:]]
                         ,facecolor='lightcoral',alpha=0.3,step='mid')

        # legend and saving figure
        ax4.set_zorder(-1)
        ax1.legend(loc="upper left", bbox_to_anchor=(1.15, 0.90)
                   ,handlelength=0,borderaxespad=0, frameon=False
                   ,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
        ax2.legend(loc="upper left", bbox_to_anchor=(1.15, 0.98)
                   ,handlelength=0,borderaxespad=0, frameon=False
                   ,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
        ax3.legend(loc="upper left", bbox_to_anchor=(1.15, 0.82)
                   ,handlelength=0,borderaxespad=0, frameon=False
                   ,fontsize=18, fancybox=False,facecolor='white',framealpha=1)
        fig.set_size_inches((17, 6))
        fig.savefig(f"{self.plotname}.png",transparent=False,
                    bbox_inches='tight')

    def _weighted_ce(self,ce_value_array,ce_error_array):
        """ calculate some weighted c/e averages

        Parameters
        ----------
        ce_value_array : list[float]
            Array of C/E values to do the maffs on
        ce_error_array : list[float]
            See above but for errors

        Returns
        -------
        weighted_ce_result : float
            weighted c/e value for requested c/e results
        weighted_ce_error : float
            weighted c/e uncert for requested c/e results
        """
        weights = []
        weighted_values = []
        for i in np.arange(len(ce_value_array)):
            if ce_value_array[i] > 0:
                weight = 1/((ce_error_array[i])**2)
                weighted_value = weight*ce_value_array[i]
                weights.append(weight)
                weighted_values.append(weighted_value)
        summed_weights = np.sum(weights)
        summed_weighted_values = np.sum(weighted_values)
        weighted_ce_result = summed_weighted_values/summed_weights
        weighted_ce_error = 1/np.sqrt(summed_weights)
        return weighted_ce_result,weighted_ce_error
    
    def run(self,calc_results,exp_results,flux_norm,flux_error,
            we_isotopes,libraries,new_order,plot_splitter,
            y_axis,we_library):
        """ read the c_results data and e_results data, 
        calculate C/E results and uncerts and plot 
        and do weighted average analysis

        Parameters
        ----------
        calc_results : str
            String name of the calculated results (without .json)
        exp_results : str
            String name of the experimental results (without .json)
        flux_norm : float
            Normalisation factor for the C/E results i.e. if you didn't
            normalise the flux properly in your fispact calculation 
            and need to do a bit more
        flux_error : float
            Fractional uncertainty on the flux estimation
        we_isotopes : list[int]
            specify isotopes you want to do weighted average analysis on 
            like [4,19] for isotopes 4-->19
        libraries : list[str]
            List of the three library names that you are using
        new_order : list[int]
            New order of the original isotopes to plot in 
            i.e. [4,3,2,1,0] will reverse an original list of 5 isotopes
        plot_splitter : int
            How many isotopes in the new_order list would you like to be
            split on the LHS of the plot
        y_axis : list[float]
            set y axis limits i.e. [0,2] for C/E = 0-->2
        we_library : str
            Name of the library to do the weighted ave analysis on 
            (must be one of 'libraries')
        """
        with open(f"{calc_results}.json") as model_results_path:
            model_results = json.load(model_results_path)

            # extract foil data and isotopes from c_results
            isotope_list = model_results.keys()
            foil_weight = [model_results[i]["foil_weight"]
                           for i in isotope_list]
            ssf = [model_results[i]["self_shielding_factor"]
                   for i in isotope_list]
            isotope_list_mathmode = [model_results[i]["mathmode_name"]
                                     for i in isotope_list]

            # get isotopic spectrum + flux uncertainties from c_results 
            spectrum_frac_u = [model_results[i]
                               ["spectrum_fractional_uncertainty"]
                               for i in isotope_list]
            spectrum_flux_frac_u = [np.sqrt(flux_error**2 + u**2)
                                    for u in spectrum_frac_u]

            # extract calculated activities from c_results 
            calc_a_1 = [model_results[i]["activities"][0]
                        for i in isotope_list]
            calc_u_1 = [model_results[i]["fractional_uncertainties"][0]
                        for i in isotope_list]
            calc_a_2 = [model_results[i]["activities"][1]
                        for i in isotope_list]
            calc_u_2 = [model_results[i]["fractional_uncertainties"][1]
                        for i in isotope_list]
            calc_a_3 = [model_results[i]["activities"][2]
                        for i in isotope_list]
            calc_u_3 = [model_results[i]["fractional_uncertainties"][2]
                        for i in isotope_list]

        with open(f"{exp_results}.json" ) as exp_results_path:
            exp_results_data = json.load(exp_results_path)

            # extract experimental activities from e_results
            # (using average activities of included peaks)
            exp_a = [np.mean(exp_results_data[key]["activities"])
                     for key in isotope_list]
            exp_u = [np.mean(exp_results_data[key]["activity_uncertainties"])
                     for key in isotope_list]

        # perform C/E calculations for the foils
        ce_results_1  = [
            (flux_norm)*
            self._c_over_e(calc_a_1,exp_a,isotope_list,foil_weight,ssf)
            [i] for i in new_order]
        ce_results_2  = [
            (flux_norm)*
            self._c_over_e(calc_a_2,exp_a,isotope_list,foil_weight,ssf)
            [i]  for i in new_order]
        ce_results_3 = [
            (flux_norm)*
            self._c_over_e(calc_a_3 ,exp_a,isotope_list,foil_weight,ssf)
            [i]  for i in new_order]
        ce_errors_1 =   [
            (flux_norm)*
            self._c_over_e(calc_a_1,exp_a,isotope_list,foil_weight,ssf)[i]
            *self._c_over_e_uncerts(calc_u_1,exp_a,exp_u,isotope_list,
                                    spectrum_flux_frac_u
                                    )[i] for i in new_order]
        ce_errors_2 =   [
            (flux_norm)*
            self._c_over_e(calc_a_2,exp_a,isotope_list,foil_weight,ssf)[i]
            *self._c_over_e_uncerts(calc_u_2,exp_a,exp_u,isotope_list,
                                    spectrum_flux_frac_u
                                    )[i] for i in new_order]
        ce_errors_3 =  [
            (flux_norm)*
            self._c_over_e(calc_a_3,exp_a,isotope_list,foil_weight,ssf)[i]
            *self._c_over_e_uncerts(calc_u_3,exp_a,exp_u,isotope_list,
                                    spectrum_flux_frac_u)[i]for i in new_order]
        spectrum_flux_frac_u = [spectrum_flux_frac_u[i] for i in new_order]
        
        #print some results
        new_isotope_list = [isotope_list_mathmode[i] for i in new_order]
        for i in range(len(new_isotope_list)):
            print(f'********* {new_isotope_list[i]} C/E results')
            print(f"{libraries[0]} value is "
                  f"{ce_results_1[i]} +- {ce_errors_1[i] }")
            print(f"{libraries[1]} value is "
                  f"{ce_results_2[i]} +- {ce_errors_2[i] }")
            print(f"{libraries[2]} value is "
                  f"{ce_results_3[i]} +- {ce_errors_3[i] }")
        
        # plot results
        self._plotter(new_order,new_isotope_list,ce_results_1,
                      ce_errors_1,ce_results_2,ce_errors_2,
                      ce_results_3,ce_errors_3,
                      spectrum_flux_frac_u,libraries,plot_splitter,
                      y_axis)
        
        # find out which library is being used for WE calculations
        we_library_index = libraries.index(we_library)
        if we_library_index == 0:
            we_ce_results = ce_results_1
            we_ce_uncerts = ce_errors_1
        if we_library_index == 1:
            we_ce_results = ce_results_2
            we_ce_uncerts = ce_errors_2
        if we_library_index == 2:
            we_ce_results = ce_results_3
            we_ce_uncerts = ce_errors_3

        # do weighted ave calcs and print 
        weighted_ave_value = self._weighted_ce(
            we_ce_results[we_isotopes[0]:we_isotopes[1]],
            we_ce_uncerts[we_isotopes[0]:we_isotopes[1]])[0]
        weighted_ave_uncert = self._weighted_ce(
            we_ce_results[we_isotopes[0]:we_isotopes[1]],
            we_ce_uncerts[we_isotopes[0]:we_isotopes[1]])[1]
        print(f"weighted {libraries[we_library_index]} C/E for "
              f"{new_isotope_list[we_isotopes[0]:we_isotopes[1]]} is "
              f"{weighted_ave_value} +- {weighted_ave_uncert}")