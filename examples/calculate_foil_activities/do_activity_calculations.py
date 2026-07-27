""" calculate experimental foil activities and reaction rates 
from a json of foil measurements and a json of calibration curves.
"""
from bfoils.activity import ActivityCalc
from datetime import datetime

# path to the directory containing the input data jsons
# the results file will also be saved here
experiment_dir = '.'

# filepath for the gamma spec measurements json
# should include dictionaries for each measured isotope including:
#  live_time of the measurement in s
#  datetime of the start of the measurement in datetime format
#  counts for each peak in decreasing order of intensity
#  raw uncertainty for each peak in the same order
#  coincidence summing correction factor for each peak in same order
#  radius of the foil in cm
#  calibration curve in the calibration json used for this measurement
#  thickness of the foil in cm
#  density of the foil in gcm3
#  element of the foil i.e. "fe"
measurements_file = 'data/measurement_data.json'

# filepath for the detector calibration curves json
# should include dictionaries for each calibration including:
#  distance from calibration points to detector in cm 
#  radius of the detector crystal in cm
#  parameters for the efficiency equation = 
#    exp(a0 + a1log(E) + a2log(E)^2+ a3log(E)^3)
#  average fractional uncertainty on the efficiency curve fit
calibration_file = 'data/calibrations.json'

# input datetime for the end of the irradiation
irradiation_end = datetime(2024,3,28, 18,17,00)

# load calibration and foil data
activities = ActivityCalc(measurements_file,experiment_dir,
                          irradiation_end,calibration_file)

# input the total irradiation time in seconds (for reaction rate calculation) 
irrad_time = 3600*2

# input the desired name of the json activity results file
e_results_name = 'e_results'

# path to xcom attenuation data folder
xcom_folder = '../../data/XCOM_new'

# get experimental results json with activities and reaction rates
# uncertainties stem from peak fitting and efficiency uncertainty
activities.calculate_activities(irrad_time,e_results_name,xcom_folder)
