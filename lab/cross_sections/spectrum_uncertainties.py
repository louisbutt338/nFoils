""" example for getting isotopic spectrum uncertainties
"""
import numpy as np
from nfoils.reaction import IsotopicSpectrumUncertainty

# set energy grid from another file
#ek=np.fromfile('../../data/energy_grids/group_structure_211.txt', sep=" ")
ek =np.linspace(5e5, 20e6, 51)

# set library (check if available in sandy)
# need internet for sandy to get library data 
library = '../../data/endf/tendl_21'  # endfb_80 jeff_33 tendl_21

# path to spectrum data txt from working dir
spectrum_file = 'data/spectra/deuteron_nov24/jendl_dli_nov24_211.txt'

# path to foil reactions data json from working dir
reaction_file = 'data/reactions/proton_march24/foil_data_fe.json'

# endf filenames if looking locally
# same order as reactions in reaction_file
endf_files = [
    #'n_049-In-115_4931.dat',
    #'n_066-Dy-164_6649.dat',
    #'n_079-Au-197_7925.dat',
    #'n_049-In-115_4931.dat',
    #'n_029-Cu-65_2931.dat',
    'n_026-Fe-56_2631.dat',
    #'n_013-Al-27_1325.dat',
    #'n_079-Au-197_7925.dat',
    #'n_041-Nb-93_4125.dat',
    #'n_028-Ni-58_2825.dat'
]

# how many values to cut off from the ends of the group structure 
# i.e. if there are no spectrum values, 
# or if MCNP changes to model mode rather than data
# p-li: [0,5] = 1e-8 MeV to 15 MeV
# d-li1 [15,35] = 1e-5 MeV to 20 MeV
cutoffs = None

# run to get some reaction-dependent spectrum uncertainties
uncertainties = IsotopicSpectrumUncertainty(ek, library)
uncertainties.get_isotopic_uncertainties(spectrum_file,reaction_file,
                                         endf_files,cutoffs)
