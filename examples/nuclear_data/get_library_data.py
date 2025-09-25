""" example file for getting foil response functions and uncertainties out
of nuclear data libraries
"""
import numpy as np
from nfoils.reaction import PostprocessReactions

# set energy grids with sandy (need sandy import)
# or make your own (need numpy import)
#ek=sandy.energy_grids.SCALE238
ek=np.fromfile('../../data/energy_grids/group_structure_175.txt',sep=" ")

# set library (check if availble in sandy)
# need internet to get library data 
library = 'jeff_33'  # endfb_71 endfb_80 jendl_40u jeff_33  tendl_21

# path to foil data json from working dir
data_file_name = 'reactions/foil_data'

# list of reaction labels for plotting (must match data in input file)
reaction_labels = [r"${}^{115}$In(n,n')", 
                   r'${}^{65}$Cu(n,p)',
                   r'${}^{56}$Fe(n,p)']

# get some response functions and nuclear data uncertainties
reactions = PostprocessReactions(ek,library)
#reactions.run_rf(data_file_name,reaction_labels)
reactions.run_stdev(data_file_name,reaction_labels)
