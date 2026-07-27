""" convert time binned rate data into a unit of your choice
using a conversion factor, then dump into a json
"""
import json
from bfoils.history import RateConversion

# set path to the rate data txt file
# should include:
#  column 1: timespan of measurement in s
#  column 2: measured count rate in /s
filename = "data/rate.txt"

# load initial rate data
history = RateConversion(filename)

# enter conversion factor to multiply by the rate values
conversion_factor = 1.7e7

# get converted data
timings,converted_data = history.get_new_data(conversion_factor)

# dump converted rates into a json
new_rates_dict = {'timings_s': [i for i in timings],
                 'new_rates':[i for i in converted_data]}
with open('new_rates.json','w',encoding='utf-8'
          ) as new_rates_json:
    json.dump(new_rates_dict,new_rates_json,
              ensure_ascii=False,indent=4)
