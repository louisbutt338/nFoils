""" example for converting time bin rate data into currents  
using conversion factor
"""
import json
import numpy as np
from bfoils.history import RateConversion

# path to the rate data txt file with
# time(s) and count rate (1/s) in the two columns
# acquisition start time was 29/11/24 13:46:35
filename = "deuteron_nov24/diamond_rate"

# initialise
history = RateConversion(filename)

# enter conversion factor for going from:
# rate data cps --> charged particle cup current uA
# here, example based on 12uA=190k counts
cps_to_current_ua = (12/190000)

# get currents on faraday cup and input second by second timing list
current_list = history.get_new_data(cps_to_current_ua)[1]
timing_list = np.ones(len(current_list))

# dump the cup currents and timings for the experiments
cup_data_dict = {'cup_timings_s': [i for i in timing_list],
                 'cup_currents_ua':[i for i in current_list]}
with open('deuteron_nov24/cup_currents2.json','w',encoding='utf-8'
          ) as cup_json:
    json.dump(cup_data_dict,cup_json,
              ensure_ascii=False,indent=4)
