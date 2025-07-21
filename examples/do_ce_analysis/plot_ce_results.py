from nfoils.ce_plotter import CEPlotter

# select directory with two calculated and experimental results files inside
experiment_directory = 'proton_march24'
#working_directory = f'/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/analysis/{experiment}'
calculated_results_file = 'c_results_unfolded_njoy_xs'
experimental_results_file = 'e_results'

# name for the C/E plot
plotname = 'root_unfolded_180725'

# set y axis size
y_axis_upper = 3.5
y_axis_lower = 0

# FLUX NORMALISATION for C results
# calculated using be7 calcs from the xs_calculator program
# FNM: mcnp/be7/p-li= 0.7(0.026985 rescale) or 0.801787(0.0309 rescale), 0.8733 for initial unfolded175/be7/p-li fispact sims 
flux_norm_mean = 1
# FPE: 0.12783 for be7/p-li 
flux_percentage_error = 0.12783 

#first weighted ave in list and last weighted ave in list for the WE calculations
first_we = 5
last_we = 14

#run 
ce_plotter = CEPlotter(experiment_directory,calculated_results_file,
                           experimental_results_file,plotname,
                           y_axis_upper,y_axis_lower,flux_norm_mean,
                           flux_percentage_error,first_we,last_we)
ce_plotter.run()
