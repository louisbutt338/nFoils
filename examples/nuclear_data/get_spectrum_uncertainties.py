""" get the average uncertainty from an input neutron spectrum over
a specified cross section range of a reaction
"""
import numpy as np
from nfoils.reaction import IsotopicSpectrumUncertainty

# set predefined sandy energy group structure
# make sure to import sandy.energy_grids if doing this 
# see sandy repo for details: https://github.com/luca-fiorito-11/sandy
#ek=sandy.energy_grids.SCALE238

# or define your own group structure and import in e.g.
ek=np.fromfile('../../data/energy_grids/group_structure_175.txt',sep=" ")

# specify nuclear data library (check the list available in sandy)
# make sure you have internet to get library data 
library = 'jeff_33'  # endfb_71 endfb_80 jendl_40u jeff_33  tendl_21

# set path from the working dir to the json with all the foil reactions
# in it with the following info:
#  MAT ENDF6 material number of parent isotope e.g. 491150
#  MT ENDF6 reaction number e.g. 102
#  density of the foil in g/cm3
#  mass of the foil in g
#  atomic mass of the foil in amu
#  abundance of the parent isotope in the foil
data_file_name = 'reactions/foil_data'

# set path from working dir to the spectrum data txt file, 
# which should include:
#  the neutron spectrum in column 1
#  and associated raw uncertainty in column 2
spectrum_file = 'spectra/example'

# input no. of values to cut off from the ends of the group structure
# i.e. if there are no spectrum values at the fast end of the group structure 
cutoff = [15,35]

# get some reaction-dependent spectrum uncertainties 
uncertainties = IsotopicSpectrumUncertainty(ek,library)
uncertainties.get_isotopic_uncertainties(spectrum_file,data_file_name,cutoff)