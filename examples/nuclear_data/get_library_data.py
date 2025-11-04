""" use NJOY exe to get reaction response functions and uncertainties out
of a specified nuclear data library. note libraries may be missing reactions 
and/or covariances
"""
import numpy as np
from nfoils.reaction import PostprocessReactions

# set energy grids with sandy (requires sandy import)
# or import your own like "group_structure_175"
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
