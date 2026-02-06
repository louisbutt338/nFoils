""" get the average uncertainty from an input neutron spectrum uncertainty 
array, for a specific set of nuclear reactions

These can be input into the c_results file when doing c/e analysis
"""
import numpy as np
from nfoils.reaction import IsotopicSpectrumUncertainty

# set predefined sandy energy group structure
# make sure to import sandy.energy_grids if doing this 
# see sandy repo for details: https://github.com/luca-fiorito-11/sandy
#ek=sandy.energy_grids.SCALE238

# or define your own group structure and import in e.g.
ek=np.fromfile('data/group_structure.txt',sep=" ")

# specify nuclear data library
# check which folders are available in data/endf
library = 'tendl_21'

# set path to a json with foil reaction data
# should include dictionaries for each reaction including:
#  MAT ENDF6 material number of parent isotope e.g. 491150
#  MT ENDF6 reaction number e.g. 102
#  density of the foil in g/cm3
#  mass of the foil in g
#  atomic mass of the foil in amu
#  abundance of the parent isotope in the foil
data_file_name = 'data/reactions'

# set path to a spectrum data txt file, 
# should include two columns:
#  column 1: the neutron spectrum in arbitrary units
#  column 2: associated raw uncertainty in the same units
spectrum_file = 'data/spectrum'

# endf files, should be saved to data/endf/{library}
# put in same order as reactions in reaction_file
endf_files = [
    'n_079-Au-197_7925.dat',
    'n_041-Nb-93_4125.dat',
    'n_028-Ni-58_2825.dat'
]

# load group structure and set library
uncertainties = IsotopicSpectrumUncertainty(ek,library)

# get the reaction-specific spectrum uncertainties 
uncertainties.get_isotopic_uncertainties(spectrum_file,data_file_name,
                                         endf_files)