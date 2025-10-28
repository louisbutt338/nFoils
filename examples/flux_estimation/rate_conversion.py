""" example for converting time bin rate data into another unit
using conversion facto
"""

import json
from nfoils.history import RateConversion

# path to the rate data txt file with
# timespan(s) and count rate (1/s) in the two columns
filename = "rate"

# initialise
history = RateConversion(filename)

#enter conversion factor to go from rate to another value
conversion_factor = 1.7e7

# get data
timings,converted_data = history.get_new_data(conversion_factor)

# dump data in a json
new_data_dict = {'timings_s': [i for i in timings],
                 'new_units_data':[i for i in converted_data]}
with open('new_data.json','w',encoding='utf-8'
          ) as new_data_json:
    json.dump(new_data_dict,new_data_json,
              ensure_ascii=False,indent=4)
