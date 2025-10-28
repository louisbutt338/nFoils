""" 
module for doing various fispact-y postprocessings
"""

import json
import actigamma as ag 
import os
import matplotlib.pyplot as plt 
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=22)


class JsonRetriever:
    """ retrieve activities at a specified time interval from fispact json. 
    notes: expand to include metastable isotopes would be useful
    """

    def __init__(self,filepath):
        """ Initialise JsonRetriever class

        Attributes
        ----------
        filepath : str
            path to the fispact json file to be analysed
        """
        # set attributes
        self.filepath = filepath

    def get_acts(self,time_interval,zai):
        """ get activities for all requested isotopes

        Parameters
        ----------
        time_interval : int
            time interval in fispact json
        zai : int
            zai number for the isotope requested

        Returns
        ----------
        isotope_activity : float
            activity of the isotope in bq
        """
        interval = time_interval-1
        try:
            with open(f'{self.filepath}.json') as json_path:
                json_file_data = json.load(json_path)
                nuclides_list = (json_file_data['inventory_data']
                                 [interval]
                                 ['nuclides'])
                nuclides_dict = {item['zai']:item for item in nuclides_list}
                if zai in nuclides_dict.keys():
                    isotope_activity = nuclides_dict[zai]['activity']
                    return isotope_activity
                else:
                    return 'isotope not calculated'
        except FileNotFoundError:
            return 'file not found'


class GammaSpectrumModel:
    """ simple gamma spectrum modeller with actigamma
    """

    def __init__(self):
        """ initialise GammaSpectrumModel class
        """

    def plot_gamma(self,isotope_name,activity):
        """for input isotope and activity:
        prints halflife + gamma energies + intensities
        then plots simple model spectrum

        Parameters
        ----------
        isotope_name : str
            name of the isotope
        activity : float
            activity of the isotope
        """
        SPECTYPE = "gamma"
        db = ag.Decay2012Database()

        print(f'{isotope_name} actigamma data:')
        # get halflife of isotope
        print('- halflife: ', f'{db.gethalflife(isotope_name)} s')
        # get gamma data of isotope
        gamma_energies = db.getenergies(isotope_name, spectype=SPECTYPE)
        print(f'- gamma energies (keV): {gamma_energies/1e3}')
        gamma_intensities = db.getintensities(isotope_name,spectype=SPECTYPE)
        print(f'- gamma intensities : {gamma_intensities}')

        # define an energy grid between 0 and 4 MeV with 5,000 bins
        grid = ag.EnergyGrid(bounds=ag.linspace(0, 4e6, 10000))
        # bin the lines appropriately using single type aggregator
        lc = ag.LineAggregator(db, grid)
        inv = ag.UnstablesInventory(data=[
            (db.getzai(isotope_name), activity),])

        hist, bin_edges = lc(inv, spectype=SPECTYPE)
        plt.figure(figsize=(12,6))
        plt.xlabel('Energy (eV)') 
        plt.ylabel('Gamma activity')
        plt.tick_params(axis='y')
        plt.xlim(0,3e6)
        plt.yscale("log")
        plt.stairs(hist,bin_edges)
        plt.savefig(os.path.join(f'{isotope_name}_gammaspec.png')
                    ,transparent=False, bbox_inches='tight')