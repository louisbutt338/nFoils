""" fit hpge efficiency functions to calibration measurements 
 with a simple monte carlo simulation
"""
from nfoils.fitting import CurveFitter

# set path to efficiency data json file
# should include dictionaries for each peak energy (keV), including:
#  efficiency
#  raw uncertainty on the efficiency
#  fractional coincidence summing correction
input_data = "data/efficiencies"

# interpolation energy range in keV, as [start,end].
# for interpolating the model fit on 
# and finding the average uncertainty over this range
interp_range = [100,1800]

# initialise class
curve_fitter = CurveFitter(input_data,interp_range)

# number of monte carlo samples you want to run
no_of_monte_carlo_samples = 100

# run the monte carlo fit and plot curve and residuals
curve_fitter.run_mc(no_of_monte_carlo_samples)