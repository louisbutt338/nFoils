""" module for using sandy/njoy for nuclear data extraction
and various postprocessings
"""

import sandy
import json
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},weight='normal',size=20)
import numpy as np
import csv

class NuclearData:
    """ class for extracting cross sections and uncertainties with sandy/njoy
    """
    def __init__(self, ek,library):
        """ Initialise class
        """

        #set parameters
        self.ek = ek
        self.library = library

    def _get_endf_file(self,material):
        """ run the sandy get_endf routine

        Parameters
        ----------
        material : float
            MAT number
        """
        endf_file = sandy.get_endf6_file(self.library, "xs", material)
        return endf_file

    def _get_errorr_data(self,material,mt_value):
        """ get covariance, standard deviation from a material endf6 file

        Parameters
        ----------
        material : float
            MAT number
        mt_value : float
            MT number
        """
        mt = [mt_value]
        try:
            # use the get_errorr function to grab cov,std data
            endf_file = self._get_endf_file(material)
            ekws = dict(ek=self.ek)
            err = endf_file.get_errorr(temperature=0,err=1,chi=False, 
                                       nubar=False,prod=False,mubar=False, 
                                       errorr_kws=ekws,
                                       verbose=False)["errorr33"]
        except:
            print(f'-----> reactions not found for MAT {material}'
                  f': MT {mt_value}')
            return [],np.array([])
        covariance = err.get_cov()
        std = covariance.get_std().reset_index().query("MT in @mt")
        std["MT"] = std["MT"].astype("category")
        std["STD"] *= 100
        stdev_array = np.array(std["STD"])
        print(f'-----> found reactions for MAT {material}: MT {mt_value}')
        return covariance,stdev_array

    def _get_gendf_data(self,material,mt_value):
        """ get groupwise cross sections from a material endf6 file

        Parameters
        ----------
        material : float
            MAT number
        mt_value : float
            MT number
        """
        mt = [mt_value]
        try:
            # use the get_gendf function to grab xs data
            endf_file = self._get_endf_file(material)
            ekws = dict(ek=self.ek,nuclide_production=True,iwt=3)
            gendf = endf_file.get_gendf(minimal_processing=True, 
                                        err=0.005, temperature=0,
                                        groupr_kws=ekws)
        except:
            print(f'-----> reactions not found for MAT {material}: MT {mt_value}')
            return [],np.array([])
        xs = gendf.get_xs(mt=mt).data.to_numpy()
        xs_array = np.array([j for i in xs for j in i])
        print(f'-----> found reactions for MAT {material}: MT {mt_value}')
        return xs_array

    def _extract_array_data(self,data_array):
        """ extract data from the top function 
        and return the stdev or xs array split by individual MT reactions

        Parameters
        ----------
        data_array : array
            Big array of multiple reactions
        """
        if data_array.size > 0:
            number_of_arrays = len(data_array)/(len(self.ek)-1)
            if number_of_arrays > 1:
                #data_array_split = np.array_split(data_array,len(data_array)/(len(self.ek)-1))
                data_array_split = np.split(data_array, number_of_arrays)
            if number_of_arrays == 1:
                data_array_split = [data_array]
            #data_array_transposed = data_array.ravel()[None]
            return data_array_split
    
    def _calculate_response_function(
            self,cross_section,density, mass,abundance,atomic_mass,thickness):
        """ calculate response function from the cross section

        Parameters
        ----------
        cross_section : array
            cross section data for the reaction
        density : float
            density of foil material in g/cm3
        mass : float
            mass of foil in g
        abundance : float
            abundance ratio of the isotope which caused the reaction
        atomic_mass : float
            atomic mass of the foil
        thickness : float
            thickness of the foil in cm
        
        """
        foil_volume = mass/density
        atom_density = (abundance * density * 1e-24 * 6.022e23)/atomic_mass
        ss_correction_factor =  1 
        #ss_correction_factor = ( (1-np.exp(-atom_density*cross_section*thickness))/ (atom_density*cross_section*thickness))
        response_function = atom_density*foil_volume*cross_section *ss_correction_factor
        response_function[np.isnan(response_function)] = 0
        return response_function

    def _unpack_datafile(self,datafile):
        """ unpack the data in the json into some lists

        Parameters
        ----------
        datafile : str
            name of json data file
        labels : list
            list of reaction labels
        """
        json_file_data = json.load(open(f'{datafile}.json'))
        material_list =   [x['mat_number'] for x in json_file_data.values()]
        mt_list =         [x['mt_value'] for x in json_file_data.values() ]
        density_list =    [x['density_gcm3'] for x in json_file_data.values()]
        mass_list =       [x['mass_g'] for x in json_file_data.values()]
        abundance_list =  [x['isotope_abundance'] for x in json_file_data.values()]
        atomic_mass_list= [x['foil_atomic_mass'] for x in json_file_data.values()]
        thickness_list =  [x['thickness_cm'] for x in json_file_data.values()]
        return (material_list,mt_list,density_list,
            mass_list,abundance_list,atomic_mass_list,
            thickness_list)
    
class PostprocessReactions(NuclearData):
    """ class for exporting and plotting response fn and uncertainty
    """
    def __init__(self, ek,library):
        """ Initialise class
        """
        super().__init__(ek,library)

    def _export_and_plot_stdev(
            self,material_list,mt_values_list,reaction_labels):
        """ export stdev data to one csv and plot uncertainty percentages

        Parameters
        ----------
        material_list : list
            list of materials
        mt_values_list : list
            list of mt values
        reaction_labels : list
            list of reaction labels for plot
        """
        # inital figure and csv setting
        open('uncertainty.csv','w').close()
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(18,8),
                                      gridspec_kw={'width_ratios': [2, 3.5]},
                                      tight_layout=True)
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(material_list))))

        # loop through specified materials and MT values
        for material,mt,reaction in zip(material_list,mt_values_list,reaction_labels):
            nuclear_data = self._get_errorr_data(material,mt)
            array_of_arrays = self._extract_array_data(nuclear_data[1])

            # plot uncertainty data
            if array_of_arrays != None:
                ek_mev = [(i/1e6) for i in self.ek]
                for mt_iterator in range(len(array_of_arrays)):
                    c=next(color)
                    ax1.stairs(array_of_arrays[mt_iterator], ek_mev,label=f'{reaction}',color=c,lw=1.5)
                    ax2.stairs(array_of_arrays[mt_iterator], ek_mev,label=f'{reaction}',color=c,lw=1.5)
                # export data to one csv file
                for xs_stdev in array_of_arrays:
                    with open('uncertainty.csv','a',newline='') as f:
                        writer=csv.writer(f,delimiter=',' )
                        writer.writerow(xs_stdev*(1/100))
            else:
                continue

        # final plotting params
        ax1.set_xlim(1e-8, 1e0)
        ax1.set_ylim(1e0,2e2)
        ax1.set_xscale('log')
        ax1.grid()
        ax2.set_xlim(1e0,18)
        ax2.set_ylim(1e0,2e2)
        ax2.tick_params(axis='y',left=False,labelleft=False)
        ax2.grid()  
        ax2.legend( loc="upper right", frameon=True,fontsize=18,fancybox=False,
                   facecolor='white',framealpha=1,ncol=3)
        fig.supylabel(r"Standard deviation ($\%$)",y=0.55)
        fig.supxlabel("Neutron energy (MeV)",y=0.03)
        fig.savefig('percentage_uncert.png')

    def _export_and_plot_rf(self,
            material_list,mt_list,density_list,mass_list,abundance_list,
            atomic_mass_list,thickness_list,labels_list):
        """ export response function data to one csv and plot 

        Parameters
        ----------
        material_list : list
            list of materials
        mt_list : list
            list of mt values
        density_list : list
            list of densities in g/cm3
        mass_list : list
            list of masses in g
        abundance_list : list
            list of isotope abundances in ratio form
        atomic_mass_list : list
            list of atomic masses
        labels_list : list
            list of reaction labels for plot
        thickness_list : list
            list of foil thicknesses in cm
        """
        # inital figure and csv setting
        open('response_function.csv','w').close()
        np.savetxt("group_structure.csv", self.ek.ravel()[None],delimiter=',')
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(18,8),
                                    gridspec_kw={'width_ratios': [2, 3.5]},
                                    tight_layout=True)
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(material_list))))

        # loop through specified materials and MT values
        for (material,mt,density,mass,abundance,
             atomic_mass,thickness,reaction_label) in zip(
                 material_list,mt_list,density_list,
                 mass_list,abundance_list,atomic_mass_list,
                 thickness_list,labels_list):
            nuclear_data = self._get_gendf_data(material,mt)
            array_of_arrays = self._extract_array_data(nuclear_data)

            # plot data
            if array_of_arrays != None:
                ek_mev = [(i/1e6) for i in self.ek]
                for m in range(len(array_of_arrays)):
                    c=next(color)
                    cross_section = array_of_arrays[m]
                    response_function = self._calculate_response_function(
                        cross_section,density, mass,abundance,
                        atomic_mass,thickness)
                    ax1.stairs(response_function, ek_mev,label=f'{reaction_label}',color=c,lw=1.5)
                    ax2.stairs(response_function, ek_mev,label=f'{reaction_label}',color=c,lw=1.5)

                # export data to one csv file
                with open('response_function.csv','a',newline='') as f:
                    writer=csv.writer(f,delimiter=',' )
                    writer.writerow(response_function)
            else:
                continue
        
        # final plotting params
        ax1.set_xlim(1e-8, 1e0)
        ax1.set_ylim(1e-11,1e3)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.grid()
        ax2.set_xlim(1e0,18)
        ax2.set_ylim(1e-11,1e3)
        ax2.set_yscale('log')
        ax2.tick_params(axis='y',left=False,labelleft=False)
        ax2.grid()  
        ax2.legend( loc="upper right", frameon=True,fontsize=18,fancybox=False,
                   facecolor='white',framealpha=1,ncol=3)
        fig.supylabel("Response function Rn(E) (cm$^2$)",y=0.55)
        fig.supxlabel("Neutron energy (MeV)",y=0.03)
        fig.savefig('response_function.png')

    def run_rf(self,datafile,labels):
        """ run for response functions

        Parameters
        ----------
        datafile : str
            name of json data file
        labels : list
            list of reaction labels
        """
        data_lists = self._unpack_datafile(datafile)
        self._export_and_plot_rf(*data_lists,labels)

    def run_stdev(self,datafile,labels):
        """ run for stdev

        Parameters
        ----------
        datafile : str
            name of json data file
        labels : list
            list of reaction labels
        """
        data_lists = self._unpack_datafile(datafile)
        self._export_and_plot_stdev(data_lists[0],data_lists[1],labels)


class IsotopicSpectrumUncertainty(NuclearData):
    """ class for getting isotopic spectrum uncertainty
    """
    def __init__(self, ek,library):
        """ Initialise class
        """
        super().__init__(ek,library)

    def _read_spectrum_uncert(self,spectrum_file):
        """ Read fractional spectrum uncertainty from txt

        Parameters
        ----------
        spectrum_file : str
            name of spectrum file
        Returns
        ----------
        frac_uncert : array
            full spectrum fractional uncertainty array
        """
        spectrum_data=np.fromfile(f'{spectrum_file}.txt',sep=" ")
        frac_uncert = np.divide( spectrum_data[1::2], spectrum_data[::2])

        #spectrum_json = json.load(open(f'{spectrum_file}.json'))
        #frac_uncert = np.array(spectrum_json['unc_value'])

        #frac_uncert[np.isnan(frac_uncert)] = 0
        #frac_uncert[frac_uncert == np.inf] = 0
        return frac_uncert

    def _normalise_xs(self,xs):
        """ Normalise cross section and output array

        Parameters
        ----------
        xs : array
            xs array in the same group structure
        Returns
        ----------
        normalised_xs : array
            normalised xs array
        """
        total_xs = sum(xs)
        normalised_xs = np.array([i/total_xs for i in xs])
        return normalised_xs

    def _isotopic_uncertainty(self,uncertainties,xs):
        """ multiply normalised cross section by an uncertainty array
        to calculate the fractional uncertainty specific to that reaction.
        must be same energy structure as used here

        Parameters
        ----------
        uncertainties : array
            Uncertainty array in the same group structure
        xs : array
            Normalised xs array in the same group structure
        """
        percent_uncertainty = np.dot(uncertainties, xs)
        return percent_uncertainty
    
    def get_isotopic_uncertainties(self,spectrum_file,datafile,cutoff):
        """ Get all the uncertainties for your reactions
        """
        data_lists = self._unpack_datafile(datafile)
        #self.ek = self.ek[:-cutoff]

        # loop through specified materials and MT values
        uncertainties = []
        for (material,mt) in zip(data_lists[0],data_lists[1]):
            nuclear_data = self._get_gendf_data(material,mt)
            array_of_arrays = self._extract_array_data(nuclear_data)

            # get the uncertainties for each
            if array_of_arrays != None:
                for m in range(len(array_of_arrays)):
                    cross_section = array_of_arrays[m]
                    norm_cross_section = self._normalise_xs(cross_section)[:-cutoff]
                    spectrum_uncert_array = self._read_spectrum_uncert(spectrum_file)[:-cutoff]
                    uncertainty = self._isotopic_uncertainty(
                        spectrum_uncert_array,norm_cross_section)
                    uncertainties.append(uncertainty)

        print('List of isotopic spectrum uncertainties:',np.array(uncertainties))
