from nfoils.activity import ActivityCalc
from datetime import datetime

# generates an e_results file to plot against some example c_results

# type cd in filename for the input data (without .json)
data_file_name = 'interspec_data'

# input datetime for the end of the irradiation
irradiation_end = datetime(2024,11,29, 14,00,00)

# workin dir for input data and to save results e..g 'example'
experiment_dir = "example"

# input instance params
get_activities = ActivityCalc(data_file_name,experiment_dir,
                              irradiation_end)

# choose which isotopes to analyse from the data json:
# 'all' or specific isotope i.e. 'Mn56'
which_isotopes = 1

# input ave fractional uncertainty on your eff curve fit 
# across peak energy measurement range
efficiency_uncert_frac = 0.05

# input the total irradiation time in s(for the reaction rate calculation) 
irrad_time = 600

# run get_activities
get_activities.run(which_isotopes,efficiency_uncert_frac,irrad_time)
