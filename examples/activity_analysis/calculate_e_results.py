""" example for calculating experimental foil activities - 
generates an e_results file to be analysed against some example c_results
using do_ce_analysis.py
"""

from nfoils.activity import ActivityCalc
from datetime import datetime

# path to dir containing input data jsons
# results are also saved here
experiment_dir = "example"

# filename for the activation data json file
data_file_name = 'interspec_data'

# filename for the detector calibration data json file
#   should include dictionaries for each calibration including:
#   distance from calibration points to detector in cm 
#   radius of the detector crystal in cm
#   parameters for the efficiency equation = 
#     exp(a0 + a1log(E) + a2log(E)^2+ a3log(E)^3)
#   average fractional uncertainty on the efficiency curve fit
calibration_file_name = 'calibration_data'

# input datetime for the end of the irradiation
irradiation_end = datetime(2024,11,29, 14,00,00)

# initialise object
get_activities = ActivityCalc(data_file_name,experiment_dir,
                              irradiation_end,calibration_file_name)

# choose which isotopes to analyse from the data json:
# 'all' or specific isotope i.e. 'Mn56'
which_isotopes = 'all'

# input the total irradiation time in s (for reaction rate calculation) 
irrad_time = 600

# request name of json results file
results_name = 'e_results'

# run get_activities
get_activities.run(which_isotopes,irrad_time,results_name)
