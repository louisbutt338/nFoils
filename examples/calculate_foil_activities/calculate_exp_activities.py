from nfoils.activity_calculator import ActivityCalc
from datetime import datetime

# choose whether you want reaction rates along with the activities
reaction_rates = True

# choose whether to analyse all FOILS isotopes ('foils') 
# or TARGET isotopes ('target') 
# or a specific isotope ('isotope')
# assumes the first two isotopes in your datafile are the target, 
# and the rest are foils (change run() if not)
automation = 'foils'

# type cd in filename for the input data (without .json)
data_file_name = 'root_data'

# input the total irradiation time (for the reaction rate calculation) 
# and the endtime of irradiation
irrad_time = (20+67+41)*60 + 30
irradiation_end = datetime(2024,3,28, 18,17,32) 
#irrad_time = 20*60
#irradiation_end = datetime(2024,11,29, 15,18,45)

# input ave fractional uncertainty on your eff curve fit 
# across peak energy measurement range
efficiency_uncert_frac = 0.0375 #endcaps=0.0375 b03/38cm=0.04422

# workin dir for input data and to save results
json_path = "proton_march24"

get_activities = ActivityCalc(reaction_rates,automation,
                           data_file_name,efficiency_uncert_frac,
                           json_path,irrad_time,irradiation_end)
get_activities.run()
