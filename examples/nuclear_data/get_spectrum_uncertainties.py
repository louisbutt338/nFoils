""" get the average uncertainty from an input neutron spectrum uncertainty 
array, over a specific cross section range of a reaction
"""
import numpy as np
from nfoils.reaction import IsotopicSpectrumUncertainty

# set predefined sandy energy group structure
# make sure to import sandy.energy_grids if doing this 
# see sandy repo for details: https://github.com/luca-fiorito-11/sandy
#ek=sandy.energy_grids.SCALE238

# or define your own group structure and import in e.g.
ek=np.fromfile('data/group_structure.txt',sep=" ")

# specify nuclear data library (check the list available in sandy)
# make sure you have internet to get library data 
library = 'jeff_33'  # endfb_71 endfb_80 jendl_40u jeff_33  tendl_21

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

# input no. of values to remove from the ends of the group structure
# during the calculation
# i.e. if there are no spectrum values at the fast neutron end, 
# use [0,10] to cut off the last 10 values
cutoff = [15,35]

# get the reaction-specific spectrum uncertainties 
uncertainties = IsotopicSpectrumUncertainty(ek,library)
uncertainties.get_isotopic_uncertainties(spectrum_file,data_file_name,cutoff)