""" convert time binned rate data into a unit of your choice
using a conversion factor, then dump into a json
"""

import json
from nfoils.history import RateConversion

# set path to the rate data txt file
# should include:
#  column 1: timespan of measurement in s
#  column 2: measured count rate in /s
filename = "rate"

# initialise
history = RateConversion(filename)

# enter conversion factor to multiply by the rate values
conversion_factor = 1.7e7

# get converted data
timings,converted_data = history.get_new_data(conversion_factor)

# dump converted data into a json
new_data_dict = {'timings_s': [i for i in timings],
                 'new_units_data':[i for i in converted_data]}
with open('new_data.json','w',encoding='utf-8'
          ) as new_data_json:
    json.dump(new_data_dict,new_data_json,
              ensure_ascii=False,indent=4)
