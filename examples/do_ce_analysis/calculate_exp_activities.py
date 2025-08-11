from nfoils.activity import ActivityCalc
from datetime import datetime

# type cd in filename for the input data (without .json)
data_file_name = 'interspec_data'

# input datetime for the end of the irradiation
#irradiation_end = datetime(2024,3,28, 18,17,32) 
irradiation_end = datetime(2024,11,29, 15,18,45)

# workin dir for input data and to save results
json_path = "deuteron_nov24"

# input instance params
get_activities = ActivityCalc(data_file_name,json_path,
                              irradiation_end)


# choose whether to analyse all FOILS isotopes ('foils') 
# or TARGET isotopes ('target') 
# or a specific isotope ('isotope')
# assumes the first two isotopes in your datafile are the target, 
# and the rest are foils (change run() if not)
automation = 'Be7'

# input ave fractional uncertainty on your eff curve fit 
# across peak energy measurement range
efficiency_uncert_frac = 0.04422 #endcaps=0.0375 b03/38cm=0.04422

# input the total irradiation time (for the reaction rate calculation) 
#irrad_time = (20+67+41)*60 + 30
irrad_time = 20*60

# run get_activities
get_activities.run(automation,efficiency_uncert_frac,irrad_time)
