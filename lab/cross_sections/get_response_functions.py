""" example for getting response functions for reactions
"""
import numpy as np
from bfoils.reaction import PostprocessReactions

# set energy grid from another file
#ek=np.fromfile('../../data/energy_grids/group_structure_175.txt', sep=" ")
#ek = np.concatenate((np.array([1e4]),np.linspace(1e5, 20e6, 100)))
ek=np.linspace(1e5, 20e6, 50)

# path to local nuclear data lib
library = '../../data/endf/tendl_21' 

# path to input reaction data file
reaction_data_file = 'data/reactions/proton_march24/foil_data_fast.json'

# list of reaction labels for plotting (must match data in input file)
reaction_labels = [
#    #r'${}^{115}$In(n,$\gamma$)',
#    r'${}^{164}$Dy(n,$\gamma$)',
#    r'${}^{197}$Au(n,$\gamma$)',
    r"${}^{115}$In(n,n')", 
    r'${}^{65}$Cu(n,p) *',
    r'${}^{56}$Fe(n,p)',
    r'${}^{27}$Al(n,$\alpha$)', 
    r'${}^{197}$Au(n,2n)',
    r'${}^{93}$Nb(n,2n)',
    r'${}^{58}$Ni(n,2n) '
]
#reaction_labels = [r'${}^{56}$Fe(n,p)']
# reaction_labels = [r'${}^{115}$In(n,$\gamma$)']

# endf filenames if looking locally
# same order as reactions in reaction_file
endf_files = [
    #'n_049-In-115_4931.dat',
    #'n_066-Dy-164_6649.dat',
    #'n_079-Au-197_7925.dat',
    'n_049-In-115_4931.dat',
    'n_029-Cu-65_2931.dat',
    'n_026-Fe-56_2631.dat',
    'n_013-Al-27_1325.dat',
    'n_079-Au-197_7925.dat',
    'n_041-Nb-93_4125.dat',
    'n_028-Ni-58_2825.dat'
]
#endf_files = ['n_026-Fe-56_2631.dat']

# get some response functions and nuclear data uncertainties
reactions = PostprocessReactions(ek,library)
reactions.run_rf(reaction_data_file,reaction_labels,endf_files)
#reactions.run_stdev(reaction_data_file,reaction_labels,endf_files)
