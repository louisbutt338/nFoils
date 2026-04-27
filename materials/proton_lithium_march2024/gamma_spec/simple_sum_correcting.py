""" example for doing simple coincidence summing correction for 
multi peak calibration data
"""
from nfoils.coincidence import SimpleCorrection

# path to the json with the multi peak and single peak source
# energy and cps data in dicts
json_path = "2024/cps/eu152_am241_b03.json"

# run
correct_eu152_results = SimpleCorrection(json_path)
correct_eu152_results.run()