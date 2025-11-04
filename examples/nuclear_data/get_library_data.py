""" use NJOY and sandy to get response functions and uncertainties for 
specified reactions in a specified nuclear data library. 
note some libraries may be missing reactions and/or covariances
"""
import numpy as np
from nfoils.reaction import PostprocessReactions

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

# input list of reaction labels for plotting, as raw/f-strings 
# these must match the reactions listed in the input json file
reaction_labels = [r"${}^{115}$In(n,n')", 
                   r'${}^{65}$Cu(n,p)',
                   r'${}^{56}$Fe(n,p)']

# get some response functions and/or nuclear data uncertainties
reactions = PostprocessReactions(ek,library)
reactions.run_rf(data_file_name,reaction_labels)
reactions.run_stdev(data_file_name,reaction_labels)
