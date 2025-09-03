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
    useful although doesn't get uncertainties so check all results
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

        Parameters
        ----------
        spectrum_file : str
            name of spectrum file

        Returns
        ----------
        isotope_activity : float
            activity of the isotope in bq
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


class GammaSpectrumModel:
    """ simple gamma spectrum modeller with actigamma
    write own example to change the instance parameters for your own
    special isotope and activity
    """
    def __init__(self,isotope_name,activity):
        """ initialise
        """
        # set params
        self.isotope_name = 'Mn56'
        self.activity = 3.82e7

    def run(self):
        """ run the thing for input isotope and activity
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
            (db.getzai(self.isotope_name), self.activity),])

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
        plt.savefig(os.path.join(f'{self.isotope_name}_gammaspec.png')
                    ,transparent=False, bbox_inches='tight')