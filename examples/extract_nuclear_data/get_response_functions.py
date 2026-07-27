""" use NJOY and sandy to get response functions and uncertainties for 
specified reactions in a specified nuclear data library. 
note some libraries may be missing reactions and/or covariances.
will also output the used energy group structure - together these
can be used for unfolding
"""
import numpy as np
from bfoils.reaction import PostprocessReactions

# set predefined sandy energy group structure
# make sure to import sandy.energy_grids if using this 
# see sandy repo for details: https://github.com/luca-fiorito-11/sandy
#ek=sandy.energy_grids.SCALE238

# or define your own group structure elsewhere and import it in
ek=np.fromfile('../../data/energy_grids/group_structure_175.txt', sep=" ")

# path to local nuclear data library.
# the relevant files must be downloaded from
# https://nds.iaea.org/public/download-endf/
library = '../../data/endf/tendl_21' 

# set path to a json with foil reaction data
# should include dictionaries for each reaction including:
#  MAT ENDF6 material number of parent isotope e.g. 491150
#  MT ENDF6 reaction number e.g. 102
#  density of the foil in g/cm3
#  mass of the foil in g
#  atomic mass of the foil in amu
#  abundance of the parent isotope in the foil
data_file_name = 'data/reactions.json'

# input list of reaction labels for plotting, as r-strings/f-strings 
# these must match the order of the reactions listed in the foil data json
reaction_labels = ["${}^{197}$Au(n,2n)", 
                   "${}^{93}$Nb(n,2n)",
                   "${}^{58}$Ni(n,2n)"]

# endf nuclear data files, should be saved in library.
# list here in the same order as reactions in reaction_file
endf_files = [
    'n_079-Au-197_7925.dat',
    'n_041-Nb-93_4125.dat',
    'n_028-Ni-58_2825.dat'
]

# load group structure and set library
reactions = PostprocessReactions(ek,library)

# get response functions and/or nuclear data uncertainties
reactions.run_rf(data_file_name,reaction_labels,endf_files)
#reactions.run_stdev(data_file_name,reaction_labels,endf_files)
