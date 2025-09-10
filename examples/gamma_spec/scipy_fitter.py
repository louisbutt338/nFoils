""" example for fitting hpge efficiency functions to data
"""

from nfoils.fitting import CurveFitter

# path to efficiency data json file
input_data = "efficiency_data"

# interpolation energy range in keV [start,end]
# for interpolating the model fit on and finding the average uncertainty
interp_range = [100,1800]

# initialise class
curve_fitter = CurveFitter(input_data,interp_range)

# number of MC samples you want to try
no_of_monte_carlo_samples = 100

# run the monte carlo fit 
curve_fitter.run_mc(no_of_monte_carlo_samples)