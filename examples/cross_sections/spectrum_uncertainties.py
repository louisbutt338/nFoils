""" example file for getting spectrum uncertainties
specific to the cross section range of a reaction
"""
import numpy as np
from nfoils.reaction import IsotopicSpectrumUncertainty

# set energy grids with sandy or make your own
#ek=sandy.energy_grids.SCALE238
ek=np.fromfile('../../data/energy_grids/group_structure_175.txt',sep=" ")

# set library (check if available in sandy)
# need internet to get library data 
library = 'jeff_33'  # endfb_71 endfb_80 jendl_40u jeff_33  tendl_21

# path to foil reactions data json from working dir
data_file_name = 'reactions/foil_data'

# path to spectrum data txt from working dir
spectrum_file = 'spectra/example'

# how many values to cut off from the ends of the group structure
# i.e. if there are no spectrum values, 
# or if MCNP changes to model mode rather than data
cutoff = [15,35]

# or get some reaction-dependent spectrum uncertainties 
uncertainties = IsotopicSpectrumUncertainty(ek,library)
uncertainties.get_isotopic_uncertainties(spectrum_file,data_file_name,cutoff)