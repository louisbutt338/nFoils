"""
module for using sandy/njoy for nuclear data extraction
and various postprocessings

L Fiorito, Nuclear data uncertainty propagation to integral responses
using SANDY, Annals of Nuclear Energy 101 (359-366) 2017 
"""

import os
import json
import numpy as np
import sandy
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]},
   weight='normal', size=20)


class NuclearData:
    """ class for extracting cross sections and uncertainties with sandy/njoy
    """
    
    def __init__(self, ek, library):
        """ Initialise NuclearData class

        Attributes
        ----------
        ek : array[float]
            array of the energy binning being used
        library : str
            nuclear data library folder to extract data from
        """

        # set attributes
        self.ek = ek
        self.library = library

    def _get_endf_file(self, material, endf_path):
        """ get endf file for specified material

        Parameters
        ----------
        material : float
            MAT number
        endf_path : str
            material endf filename for locally saved file

        Returns
        -------
        endf_file : sandy endf object
            Endf file for the specified material
            See sandy documentation for more details
        """
        try:
            # get endf file from local endf files
            file_path = f'{self.library}/{endf_path}'
            endf_file = sandy.Endf6.from_file(file_path)

            # get endf file from online IAEA database
            #endf_file = sandy.get_endf6_file(self.library, "xs", material)

        except FileNotFoundError:
            print('cannot find endf file locally')

        return endf_file

    def _get_errorr_data(self, material, mt_value, endf_file):
        """ get covariance, standard deviation from a material endf6 file
        or return empty arrays (fiorito)

        Parameters
        ----------
        material : int
            MAT number
        mt_value : int
            MT number
        endf_file : str
            material endf filename for locally saved file

        Returns
        -------
        cov_array : array[float]
            big array of all covariances for the input material and MTs
        stdev_array : array[float]
            big array of all uncertainties for the input material and MTs
        """
        mt = [mt_value]
        try:
            # use the get_errorr function to grab cov,std data
            endf_file = self._get_endf_file(material,endf_file)
            ekws = dict(ek=self.ek)
            err = endf_file.get_errorr(temperature=300, err=1, chi=False,
                                       nubar=False, prod=False, mubar=False,
                                       errorr_kws=ekws,
                                       verbose=False)["errorr33"]
            covariance = err.get_cov()
            std = covariance.get_std().reset_index().query("MT in @mt")
            std["MT"] = std["MT"].astype("category")
            std["STD"] *= 100
            cov_array = np.array(covariance)
            stdev_array = np.array(std["STD"])
            print(f'-----> found reactions for MAT {material}: MT {mt}')
            return cov_array, stdev_array

        except (KeyError, ValueError):
            # if above doesn't work return empty arrays so code
            # keeps running for other reactions
            print(f'-----> reactions not found for MAT {material}'
                  f': MT {mt}')
            return np.array([]),np.array([])

    def _get_gendf_data(self, material, mt_value, endf_file):
        """ get groupwise cross sections from a material endf6 file
        or return empty array (fiorito)

        Parameters
        ----------
        material : float
            MAT number
        mt_value : float
            MT number
        endf_file : str
            material endf filename for locally saved file

        Returns
        -------
        xs_array : array[float]
            big array of cross sections for the material and MTs
        """
        mt = [mt_value]
        try:
            # use the get_gendf function to grab xs data
            endf_file = self._get_endf_file(material,endf_file)
            ekws = dict(ek=self.ek, nuclide_production=True, iwt=3)
            gendf = endf_file.get_gendf(minimal_processing=True,
                                        err=0.005, temperature=300,
                                        groupr_kws=ekws)
            xs = gendf.get_xs(mt=mt).data.to_numpy()
            xs_array = np.array([j for i in xs for j in i])
            print(f'-----> found reactions for MAT {material}: MT {mt}')
            return xs_array

        except ValueError:
            # if above doesn't work return empty arrays
            print(f'-----> reactions not found for MAT {material}:'
                  f' MT {mt}')
            return np.array([])

    def _extract_array_data(self, data_array):
        """ extract data from the top function and return the stdev
        or xs array split by individual MT reactions

        Parameters
        ----------
        data_array : array
            Big array of multiple reactions

        Returns
        -------
        data_array_split : array[array[float]]
            split the big array into an array of arrays for each reaction
        """
        if data_array.size > 0:
            number_of_arrays = len(data_array)/(len(self.ek)-1)
            if number_of_arrays > 1:
                data_array_split = np.split(data_array, number_of_arrays)
            if number_of_arrays == 1:
                data_array_split = [data_array]
            return data_array_split
    
    def _calculate_response_function(self, cross_section, density, mass,
                                     abundance, atomic_mass):
        """ calculate response function from the cross section

        Parameters
        ----------
        cross_section : array[float]
            cross section data for the reaction
        density : float
            density of foil material in g/cm3
        mass : float
            mass of foil in g
        abundance : float
            abundance ratio of the isotope which caused the reaction
        atomic_mass : float
            atomic mass of the foil

        Returns
        -------
        response_function : array[float]
            response function array for the reaction
        """
        foil_volume = mass/density
        atom_density = (abundance * density * 1e-24 * 6.022e23)/atomic_mass
        response_function = (atom_density * foil_volume * cross_section)
        response_function[np.isnan(response_function)] = 0
        return response_function

    def _unpack_datafile(self, datafile):
        """ unpack the data in the json into some lists

        Parameters
        ----------
        datafile : str
            name of json data file
        labels : list[str]
            list of reaction labels

        Returns
        -------
        material_list : list[int]
            list of ENDF material numbers
        mt_list : list[int]
            list of reaction (MT) numbers
        density_list : list[float]
            list of densities (g/cm3) for foils
        mass_list : list[float]
            list of masses (g) for foils
        abundance_list : list[float]
            list of abundances for parent isotopes
        atomic_mass_list : list[float]
            list of atomic masses for foils
        """
        with open(datafile) as json_file:
            json_file_data = json.load(json_file)
            material_list = [x['mat_number'] for x in json_file_data.values()]
            mt_list = [x['mt_value'] for x in json_file_data.values()]
            density_list = [x['density_gcm3'] for x in json_file_data.values()]
            mass_list = [x['mass_g'] for x in json_file_data.values()]
            abundance_list = [x['isotope_abundance']
                              for x in json_file_data.values()]
            atomic_mass_list = [x['foil_atomic_mass']
                                for x in json_file_data.values()]
        return (material_list, mt_list, density_list,
                mass_list, abundance_list, atomic_mass_list)


class PostprocessReactions(NuclearData):
    """ class for exporting and plotting response function and uncertainty 
    """

    def __init__(self, ek, library):
        """ Initialise PostprocessReactions class (inherits from NuclearData)

        Attributes
        ----------
        ek : array[float]
            array of the energy binning being used
        library : str
            nuclear data library folder to extract data from
        """
        super().__init__(ek, library)

        # make new folder to save the responses/uncertainties in
        self.results_folder = os.path.join(f"for_unfolding")
        os.makedirs(self.results_folder,exist_ok=True)

    def _export_and_plot_stdev(self, material_list, mt_values_list,
                               reaction_labels, endf_files):
        """ export stdev data to one txt and plot uncertainty percentages

        Parameters
        ----------
        material_list : list[int]
            list of endf material numbers
        mt_values_list : list[int]
            list of endf mt values
        reaction_labels : list[str]
            list of reaction labels for plot
        endf_files : list[str]
            Names for the endf files saved locally
        """
        # inital figure and txt list settings
        big_response_function_uncert_list = []
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8),
                                       gridspec_kw={'width_ratios': [2, 3.5]},
                                       tight_layout=True)
        #fig,ax2 = plt.subplots(1, figsize=(12, 8),tight_layout=True)
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(material_list))))

        # loop through specified materials and MT values
        for (material, mt, 
             reaction, e_file) in zip(material_list, mt_values_list,
                                        reaction_labels,endf_files):
            nuclear_data = self._get_errorr_data(material, mt, e_file)
            array_of_arrays = self._extract_array_data(nuclear_data[1])

            # plot uncertainty data
            if array_of_arrays is not None:
                ek_mev = [(i/1e6) for i in self.ek]
                for std in array_of_arrays:
                    c = next(color)
                    ax1.stairs(std,ek_mev,label=f'{reaction}',color=c,lw=1.5)
                    ax2.stairs(std,ek_mev,label=f'{reaction}',color=c,lw=1.5)
                    
                    # export data to one csv file
                    #with open('response_function_uncertainty.csv', 'a', newline='') as f:
                    #    writer = csv.writer(f, delimiter=',')
                    #    writer.writerow(std*(1/100))
                    # add data to txt file list
                    big_response_function_uncert_list.append(std)
            else:
                continue

        # export data to txt file
        np.savetxt(f"{self.results_folder}/response_matrix_uncertainties.txt",
                   big_response_function_uncert_list,
                   delimiter=',')

        # final plotting params
        ax1.set_xlim(1e-8, 1e0)
        ax1.set_ylim(1e0, 2e2)
        ax1.set_xscale('log')
        ax1.grid()
        ax2.set_xlim(1, 20)
        ax2.set_ylim(1e0, 2e2)
        ax2.tick_params(axis='y', left=False, labelleft=False)
        ax2.grid()
        ax2.legend(loc="upper right", frameon=True, fontsize=18,
                   fancybox=False, facecolor='white', framealpha=1,
                   ncol=3)
        fig.supylabel(r"Standard deviation ($\%$)", y=0.55)
        fig.supxlabel("Neutron energy (MeV)", y=0.03)
        fig.savefig(f'percentage_uncert.png')

    def _export_and_plot_rf(self, material_list, mt_list, density_list,
                            mass_list, abundance_list, atomic_mass_list,
                            labels_list, endf_files):
        """ export response function data to one txt and plots RFs

        Parameters
        ----------
        material_list : list[int]
            list of endf material numbers
        mt_list : list[int]
            list of endf mt values
        density_list : list[float]
            list of foil densities in g/cm3
        mass_list : list[float]
            list of foil masses in g
        abundance_list : list[float]
            list of parent isotope abundances in ratio form
        atomic_mass_list : list[float]
            list of atomic masses
        labels_list : list[str]
            list of reaction labels for plot
        endf_files : list[str]
            Names for the endf files saved locally
        """
        # inital figure and txt list setting
        big_response_functions_list = []
        np.savetxt(f"{self.results_folder}/group_structure.txt",
                   self.ek.ravel()[None],delimiter=',')
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8),
                                       gridspec_kw={'width_ratios': [2, 3.5]},
                                       tight_layout=True)
        #fig,ax2 = plt.subplots(1, figsize=(12, 8),tight_layout=True)
        color = iter(plt.cm.rainbow(np.linspace(0, 1, len(material_list))))

        # loop through specified materials and MT values
        for (material, mt, density, mass, abundance,
             atomic_mass, reaction_label, e_file) in zip(
                 material_list, mt_list, density_list,
                 mass_list, abundance_list, atomic_mass_list,
                 labels_list,endf_files):
            nuclear_data = self._get_gendf_data(material, mt,e_file)
            array_of_arrays = self._extract_array_data(nuclear_data)

            # plot data
            if array_of_arrays is not None:
                ek_mev = [(i/1e6) for i in self.ek]
                for xs in array_of_arrays:
                    c = next(color)
                    response_function = self._calculate_response_function(
                        xs, density, mass, abundance,
                        atomic_mass)
                    ax1.stairs(response_function, ek_mev,
                               label=f'{reaction_label}', color=c, lw=1.5)
                    ax2.stairs(response_function, ek_mev,
                               label=f'{reaction_label}', color=c, lw=1.5)

                ## export data to one csv file
                #with open('response_function.csv', 'a', newline='') as f:
                #    writer = csv.writer(f, delimiter=',')
                #    writer.writerow(response_function)
                # add data to txt file list
                big_response_functions_list.append(response_function)
            else:
                continue

        # export data to txt file
        np.savetxt(f"{self.results_folder}/response_matrix.txt",
                   big_response_functions_list,
                   delimiter=',')

        # final plotting params
        ax1.set_xlim(1e-8, 1e0)
        ax1.set_ylim(1e-11, 1e3)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.grid()
        ax2.set_xlim(1, 20)
        ax2.set_ylim(1e-11, 1e3)
        ax2.set_yscale('log')
        ax2.tick_params(axis='y', left=False, labelleft=False)
        ax2.grid()
        ax2.legend(loc="upper right", frameon=True, fontsize=18,
                   fancybox=False, facecolor='white', framealpha=1, ncol=3)
        fig.supylabel("Response function (cm$^2$)", y=0.55)
        fig.supxlabel("Neutron energy (MeV)", y=0.03)
        fig.savefig('response_function.png')

    def run_rf(self, datafile, labels,endf_files):
        """ extract response functions, dump in txt format for unfolding
        and plot

        Parameters
        ----------
        datafile : str
            name of json reaction data file
        labels : list[str]
            list of the raw string mathmode reaction labels, 
            matching the reactions in the datafile
        endf_files : list[str]
            Names for the endf files saved locally
        """
        data_lists = self._unpack_datafile(datafile)
        self._export_and_plot_rf(*data_lists, labels,endf_files)

    def run_stdev(self, datafile, labels, endf_files):
        """ extract standard deviations, dump in txt format for unfolding
        and plot

        Parameters
        ----------
        datafile : str
            name of json reaction data file
        labels : list[str]
            list of the raw string mathmode reaction labels,
            matching the reactions in the datafile
        endf_files : list[str]
            Names for the endf files saved locally
        """
        data_lists = self._unpack_datafile(datafile)
        self._export_and_plot_stdev(data_lists[0], data_lists[1],
                                    labels,endf_files)


class IsotopicSpectrumUncertainty(NuclearData):
    """ class for getting isotopic spectrum uncertainty
    """

    def __init__(self, ek, library):
        """ Initialise IsotopicSpectrumUncertainty class
        (inherits from NuclearData)

        Attributes
        ----------
        ek : array[float]
            array of the energy binning being used
        library : str
            nuclear data library folder to extract data from
        """
        super().__init__(ek, library)

    def _read_spectrum_uncert(self, spectrum_file):
        """ Read fractional spectrum uncertainty from txt file

        Parameters
        ----------
        spectrum_file : str
            name of spectrum txt file

        Returns
        ----------
        frac_uncert : array[float]
            full spectrum fractional uncertainty array
        """
        spectrum_data = np.fromfile(spectrum_file, sep=" ")
        flux_vals = spectrum_data[::2]
        uncert_vals = spectrum_data[1::2]
        frac_uncert = np.divide(uncert_vals, flux_vals)
        return frac_uncert

    def _normalise_xs(self, xs):
        """ Normalise cross section and output array

        Parameters
        ----------
        xs : array[float]
            xs array in the same group structure

        Returns
        ----------
        normalised_xs : array[float]
            normalised xs array
        """
        total_xs = sum(xs)
        normalised_xs = np.array([i/total_xs for i in xs])
        return normalised_xs

    def _isotopic_uncertainty(self, uncertainties, xs):
        """ multiply normalised cross section by an uncertainty array
        to calculate the fractional uncertainty specific to that reaction.
        must be same energy structure as used here

        Parameters
        ----------
        uncertainties : array[float]
            Uncertainty array in the same group structure
        xs : array[float]
            Normalised xs array in the same group structure

        Returns
        ----------
        percent_uncertainty : float
            percentage uncertainty for a cross section
            and a spectrum uncert array
        """
        percent_uncertainty = np.dot(uncertainties, xs)
        return percent_uncertainty
    
    def get_isotopic_uncertainties(self, spectrum_file, datafile,
                                   endf_files,cutoffs=None):
        """ Get all the uncertainties for your reactions
        and prints the results

        Parameters
        ----------
        spectrum_file : str
            Name of the spectrum txt file
        datafile : str
            Name of the json datafile
        endf_files : list[str]
            Names for the endf files saved locally
        cutoffs : list[int]
            how many vals to cut off at the end of the group structure
            if there is no spectrum there i.e. [1,5] cuts off 1 from the 
            bottom and 5 from the top of the group structure
        """
        # unpack the data
        data_lists = self._unpack_datafile(datafile)

        # do cutoffs from the energy grid if requested
        if cutoffs==None:
            spectrum_uncert_array = (self._read_spectrum_uncert
                                     (spectrum_file))
        else:
            self.ek = self.ek[cutoffs[0]:-cutoffs[1]]
            spectrum_uncert_array = (self._read_spectrum_uncert
                                     (spectrum_file)
                                     [cutoffs[0]:-cutoffs[1]])

        # loop through specified materials and MT values
        uncertainties = []
        materials = []
        mts = []
        for (material,mt,endf_file) in zip(data_lists[0],
                                           data_lists[1],endf_files):
            nuclear_data = self._get_gendf_data(material, mt,endf_file)
            array_of_arrays = self._extract_array_data(nuclear_data)

            # get the uncertainties for each
            if array_of_arrays is not None:
                for xs in array_of_arrays:
                    norm_cross_section = (self._normalise_xs(xs))
                    uncertainty = self._isotopic_uncertainty(
                        spectrum_uncert_array,norm_cross_section)
                    uncertainties.append(uncertainty)
                    materials.append(material)
                    mts.append(mt)

        for (i,j,k) in zip(materials, mts,uncertainties):
            print(f'mat={i},mt={j} fractional uncertainty is {k}')
