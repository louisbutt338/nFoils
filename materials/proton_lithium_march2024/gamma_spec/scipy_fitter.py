""" example for fitting hpge efficiency functions to experimental
efficiency measurements. measurements should include 
coincidence summing corrections if required
"""
from nfoils.fitting import CurveFitter

# path to efficiency data json file
input_data = "2024/efficiencies/b03_1cm.json"  # g11_pt4cm

# energy range in keV for fitting the data in montecarlo
# and finding the average uncertainty. [start,end]
interp_range = [100,1800]

# initialise class
curve_fitter = CurveFitter(input_data,interp_range)

# number of MC samples you want to try
no_of_monte_carlo_samples = 1000

# run the monte carlo fit 
curve_fitter.run_mc(no_of_monte_carlo_samples)