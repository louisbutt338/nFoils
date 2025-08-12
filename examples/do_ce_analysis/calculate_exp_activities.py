from nfoils.activity import ActivityCalc
from datetime import datetime

# type cd in filename for the input data (without .json)
data_file_name = 'interspec_data'

# input datetime for the end of the irradiation
irradiation_end = datetime(2024,11,29, 14,00,00)

# workin dir for input data and to save results e..g 'test'
experiment_dir = "../do_ce_analysis"

# input instance params
get_activities = ActivityCalc(data_file_name,experiment_dir,
                              irradiation_end)


# choose to analyse the first two isotopes in the data json ('target') 
# or all other isotopes ('foils')
# or a specific isotope ('isotope')
# feel free to change ActivityCalc.run for your own purpose
automation = 'Be7'

# input ave fractional uncertainty on your eff curve fit 
# across peak energy measurement range
efficiency_uncert_frac = 0.05

# input the total irradiation time in s(for the reaction rate calculation) 
irrad_time = 600

# run get_activities
get_activities.run(automation,efficiency_uncert_frac,irrad_time)
