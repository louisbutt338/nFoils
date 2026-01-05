""" do a simple coincidence summing correction for
multi peak calibration source measurements
"""
from nfoils.coincidence import SimpleCorrection

# set path to json file with multi peak and single peak measurements inside
# should include dictionaries for "multi_source_near", "multi_source_far",
# "single_source_near", "single_source_far"
# each dictionary should include:
#  "(peak energy in keV)": (peak counts per second)
json_path = "data/peak_measurements"

# get correction factors for each of the peak energies
correct_eu152_results = SimpleCorrection(json_path)
correct_eu152_results.run()