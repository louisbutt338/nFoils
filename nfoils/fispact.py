""" 
module for doing various fispact-y things
"""

import json
import actigamma as ag 
import os
import matplotlib.pyplot as plt 
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=22)


class JsonRetriever:
    """ retrieve activities at a specified time interval from fispact json
    useful although doesn't get uncerts so check all results
    """
    def __init__(self,filepath,interval,zai,isotope_state):
        """ Initialise class
        """
        # set the parameters
        self.filepath = filepath
        self.interval = interval-1
        self.zai = zai
        self.state = isotope_state

    def run(self):
        """ run for all requested isotopes
        """
        try:
            json_file_data = json.load(open(f'{self.filepath}'))
        except FileNotFoundError:
            return 'file not found'
        nuclides_list = (json_file_data['inventory_data'][self.interval]
                         ['nuclides'])
        nuclides_dict = {item['zai']:item for item in nuclides_list}
        if self.zai in nuclides_dict.keys():
            isotope_activity = nuclides_dict[self.zai]['activity']
            return isotope_activity
        else:
            return 'isotope not calculated'


class ActivitiesCollector:
    """ class to collect activities from a json
    """
    def __init__(self,library,fispact_results_folder):
        """ initialise class
        """
        # set parameters
        self.library = 'tendl21'
        self.fispact_results_folder = (
            '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/'
            f'uBB/040924_foils_fe_flux_analysis/{self.library}')
        self.materials = [
            'au', 'al', 'fe', 'in', 'nb', 'ni', 'rh',
            'sc', 'y', 'dy', 'cd','cu']

    def run(self):
        """ run for all materials requested 
        """
        for m in self.materials:
            calculated_result_json = json.load(
                open('{}/uBB_{}_cell12.json'.format(
                    self.fispact_results_folder,m)))
            activity = calculated_result_json['neutron_cell_flux'][3]['value']
            activity_uncert = (calculated_result_json['neutron_cell_flux']
                               [3]['value'])
            print(activity,activity_uncert)


class JsonPlotter:
    """ class to plot activities from a json
    """
    def __init__(self,directory,json_name):
        """ initialise
        """
        # set params
        self.directory = directory
        self.json_name =  json_name

    def _data_extraction(self,inventory_data):
        """ extract data 

        Parameters
        ----------
        inventory_data : dict
            Json dictionary of data
        """
        timestep_array = []  # in days
        total_activity_array = []  # in Bq
        total_activity_normalised_array = []  # in Bq/g
        total_dose_array = []  # in Sv/hr
        for timestep in range(0,len(inventory_data)):
            if inventory_data[timestep]['irradiation_time'] != 0:
                timestep_array.append(
                    inventory_data[timestep]['cooling_time']/(3600*24))
                total_activity_array.append(
                    inventory_data[timestep]['total_activity'])
                total_activity_normalised_array.append(
                    inventory_data[timestep]['total_activity']
                    /(1e3*inventory_data[timestep]['total_mass']))
                total_dose_array.append(
                    inventory_data[timestep]['dose_rate']['dose'])
        return (timestep_array,total_activity_array,
                total_activity_normalised_array,total_dose_array)

    def _activity_plotter(self,timestep_array,total_activity_normalised_array):
        """ plot acitivity in bq/g from the json

        Parameters
        ----------
        timstep_array : array
            Array of timesteps
        total_activity_normalised_array : array
            Array of activities
        """
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Decay time (days)') 
        ax1.set_ylabel('Activity (Bq/g)')
        ax1.tick_params(axis='y')
        ax1.set_xlim(1e-3,1e3)
        ax1.set_xscale("log")
        #ax1.set_ylim(1e0,1e6)
        ax1.set_yscale("log")
        ax1.plot(timestep_array, total_activity_normalised_array ,
                 'k-' ,     linewidth=1.5)
        ax1.axhline(y=1e5, ls='-', c='green', lw=1.5)
        ax1.axvline(x=0.04167, ls='--', c='grey')
        ax1.axvline(x=1e0, ls='--', c='grey')
        ax1.axvline(x=30, ls='--', c='grey')
        ax1.text(1.02, 0.9, 'Approx. background', 
                 transform=ax1.transAxes, fontsize=12, c='green')
        ax1.text(0.2, 1.02, '1 hour', 
                 transform=ax1.transAxes, fontsize=12, c='grey')
        ax1.text(0.44, 1.02, '1 day', 
                 transform=ax1.transAxes, fontsize=12, c='grey')
        ax1.text(0.654, 1.02, '1 month',
                 transform=ax1.transAxes, fontsize=12, c='grey')
        ax1.grid(which='major')
        fig.set_size_inches((8, 6))
        fig.savefig(os.path.join(
            self.directory, 'total_activity_{}.png'.format(self.json_name))
            ,transparent=False, bbox_inches='tight')

    def _dose_plotter(self,timestep_array,total_dose_array):
        """ plot isotopic dose over time from the json
        Parameters
        ----------
        timstep_array : array
            Array of timesteps
        total_dose_array : array
            Array of doses
        """
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Decay time (days)') 
        ax1.set_ylabel('Dose (Sv/h)')
        ax1.tick_params(axis='y')
        ax1.set_xlim(1e-3,1e3)
        ax1.set_xscale("log")
        #ax1.set_ylim(1e-20,1e-5)
        ax1.set_yscale("log")
        ax1.plot(timestep_array, total_dose_array,'k-' ,    linewidth=1.5)
        ax1.axhline(y=1e-6, ls='-', c='green', lw=1.5)
        ax1.axvline(x=0.04167, ls='--', c='grey')
        ax1.axvline(x=1e0, ls='--', c='grey')
        ax1.axvline(x=30, ls='--', c='grey')
        ax1.text(1.02, 0.9, 'Approx. background', 
                 transform=ax1.transAxes, fontsize=12, c='green')
        ax1.text(0.2, 1.02, '1 hour',
                 transform=ax1.transAxes, fontsize=12, c='grey')
        ax1.text(0.44, 1.02, '1 day',
                 transform=ax1.transAxes, fontsize=12, c='grey')
        ax1.text(0.654, 1.02, '1 month',
                 transform=ax1.transAxes, fontsize=12, c='grey')
        ax1.grid(which='major')
        fig.set_size_inches((8, 6))
        fig.savefig(os.path.join(self.directory, 
                                 'total_dose_{}.png'.format(self.json_name))
                    , transparent=False, bbox_inches='tight')

    def run(self):
        """ run both plottings
        """
        dictionary = json.load(open('{}/{}.json'.format(self.directory,
                                                        self.json_name), 'r'))
        inventory_data = dictionary['inventory_data']
        #print(inventory_data[0]['nuclides'][0])
        print('********************',
              'plotting activity in Bq/g and dose in Sv/hr',
              '********************')
        print('note: FISPACT {} dose calculation was used'.format(
            inventory_data[0]['dose_rate']['type']))
        self._activity_plotter()
        self._dose_plotter()


class GammaSpectrumModel:
    """ simple gamma spectrum modeller with actigamma
    """
    def __init__(self,folder_path,isotope_name,activity):
        """ initialise
        """
        # set params
        self.folder_path = (
            '/Users/ljb841@student.bham.ac.uk/fispact/'
            'WORKSHOP/uBB/analysis/fispact_gammaspec')
        self.isotope_name = 'Mn56'
        self.activity = 3.82e7

    def run(self):
        """ run the thing
        """
        SPECTYPE = "gamma"
        db = ag.Decay2012Database()

        # get halflife of isotope
        print(db.gethalflife(self.isotope_name))
        # get gamma lines of isotope
        print(db.getenergies(self.isotope_name, spectype=SPECTYPE))
        print(db.getintensities(self.isotope_name,spectype=SPECTYPE))

        # define an energy grid between 0 and 4 MeV with 5,000 bins
        grid = ag.EnergyGrid(bounds=ag.linspace(0, 4e6, 10000))
        # bin the lines appropriately using single type aggregator
        lc = ag.LineAggregator(db, grid)
        inv = ag.UnstablesInventory(data=[
            (db.getzai(self.isotope_name), self.activity),
        ])

        hist, bin_edges = lc(inv, spectype=SPECTYPE)
        plt.figure(figsize=(12,6))
        plt.xlabel('Energy (eV)') 
        plt.ylabel('gamma activity')
        plt.tick_params(axis='y')
        plt.xlim(0,3e6)
        #plt.set_xscale("log")
        plt.ylim(1e-5,1e15)
        plt.yscale("log")
        plt.stairs(hist,bin_edges)
        plt.savefig(os.path.join(self.folder_path, '{}_gammaspec.png'.format(
            self.isotope_name))
            ,transparent=False, bbox_inches='tight')