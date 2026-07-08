""" calculate simulated foil reaction rates 
from a json of simulation data, spectrum and response matrix.
a c_results file is generated, which can be analysed against the
e_results json file with do_ce_analysis.py
"""
from nfoils.activity import ActivitySim

# path to the directory containing the input data jsons
# the results file will also be saved here
experiment_dir = '.'

# simulation data file path
simulation_file = 'data/simulation_data'

# initialise
activities = ActivitySim(experiment_dir,simulation_file)

# choose which isotopes to calculate activity for from the data json
# 'all' or specific isotope i.e. 'Mn56'
which_isotopes = 'all'

# input the total irradiation time in seconds (for reaction rate calculation) 
irrad_time = 3600*2

# name for the output c results file
c_results_name = 'c_results_tendl_test'

# paths to the reponse matrix and associated uncertainties
response_matrix = '../nuclear_data/results/response_matrix.txt'
response_matrix_uncerts = '../nuclear_data/results/response_matrix_uncertainties.txt'

# path to the neutron spectrum in the same group structure
spectrum = '../nuclear_data/data/response_matrix.txt'

# do simulated activities with response functions from the reaction module
activities.simulate_activities(which_isotopes,irrad_time,c_results_name,
                               response_matrix,response_matrix_uncerts,spectrum) 

