""" example for doing simple coincidence summing correction for 
multi peak calibration data
"""
from nfoils.coincidence import SimpleCorrection

# path to the json with the multi peak and single peak source
# energy and cps data in dicts
json_path = "coincidence_data"

# run
correct_eu152_results = SimpleCorrection(json_path)
correct_eu152_results.run()