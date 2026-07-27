""" example for doing simple coincidence summing correction for 
multi peak calibration data
"""
from bfoils.coincidence import SimpleCorrection

# path to the json with the multi peak and single peak source
# energy and cps data in dicts
json_path = "nov24/cps/eu152_am241_b03"

# run
correct_eu152_results = SimpleCorrection(json_path)
correct_eu152_results.run()