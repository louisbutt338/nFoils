""" example for calculating foil activities
"""
from datetime import datetime
from nfoils.activity import ActivityCalc

# workin dir for input data and to save results e.g 'test'
experiment_dir = "proton_march24"

# type cd in filename for the input data 
data_file_name = 'root_data.json'

# filename for the detector calibration data file 
# should include dictionaries for each calibration including:
# distance from calibration points to detector in cm 
# radius of the detector crystal in cm
# parameters for the efficiency equation = 
#   exp(a0 + a1log(E) + a2log(E)^2+ a3log(E)^3)
# average fractional uncertainty on the efficiency curve fit
calibration_file_name = 'calibration_data.json'

# input datetime for the end of the irradiation
irradiation_end = datetime(2024,3,28, 18,17,32)

# set instance params
get_activities = ActivityCalc(data_file_name,experiment_dir,
                              irradiation_end,calibration_file_name)


# choose which isotopes to analyse from the data json using strings:
# 'all' or specific isotope i.e. 'Mn56'
# OR if you want to skip the first two isotopes in the file, put 2
# if you want to skip the last two isotopes from your file, put -2
which_isotopes = 2

# input the total irradiation time in s (for the reaction rate calculation) 
irrad_time = (20+67+41)*60 + 30

# request name of json results file
results_name = 'e_results'

# run get_activities
get_activities.calculate_activities(irrad_time,results_name,'../../../data/XCOM_new',
                                    which_isotopes)
