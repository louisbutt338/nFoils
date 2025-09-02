from nfoils.fitting import CurveFitter

#select folder and dataset
input_data_path = "../gamma_spec"
input_data_filename = "endcap_data"

# energy range in keV for fitting the data in montecarlo 
# and finding the average uncertainty
interpolation_range_start = 100
interpolation_range_end = 1800

# number of MC samples you want to try
no_of_monte_carlo_samples = 100

# run
curve_fitter = CurveFitter(input_data_path,input_data_filename,
                           interpolation_range_start,interpolation_range_end,
                           no_of_monte_carlo_samples)
curve_fitter.run_mc()