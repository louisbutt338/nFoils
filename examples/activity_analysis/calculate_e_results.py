""" calculate experimental foil activities and reaction rates 
from a json of foil measurements and a json of calibration curves.
an e_results file will be generated which can be analysed against the example
c_results json file using do_ce_analysis.py
"""

from nfoils.activity import ActivityCalc
from datetime import datetime

# path to the directory containing the input data jsons
# the results file will also be saved here
experiment_dir = "example"

# filename for the gamma spec measurements json
#   should include dictionaries for each measured isotope including:
#   live_time of the measurement in s
#   datetime of the start of the measurement in datetime format
#   counts for each peak in decreasing order of intensity
#   raw uncertainty for each peak in the same order
#   coincidence summing correction factor for each peak in same order
#   radius of the foil in cm
#   calibration curve in the calibration json used for this measurement
#   thickness of the foil in cm
#   density of the foil in gcm3
#   element of the foil i.e. "fe"
data_file_name = 'interspec_data'

# filename for the detector calibration curves json
#   should include dictionaries for each calibration including:
#   distance from calibration points to detector in cm 
#   radius of the detector crystal in cm
#   parameters for the efficiency equation = 
#     exp(a0 + a1log(E) + a2log(E)^2+ a3log(E)^3)
#   average fractional uncertainty on the efficiency curve fit
calibration_file_name = 'calibration_data'

# input datetime for the end of the irradiation
irradiation_end = datetime(2024,11,29, 14,00,00)

# initialise class
get_activities = ActivityCalc(data_file_name,experiment_dir,
                              irradiation_end,calibration_file_name)

# choose which isotopes to calculate activity for from the data json
# 'all' or specific isotope i.e. 'Mn56'
which_isotopes = 'all'

# input the total irradiation time in seconds (for reaction rate calculation) 
irrad_time = 600

# input the desired name of the json activity results file
results_name = 'e_results'

# run the code
get_activities.run(which_isotopes,irrad_time,results_name)
