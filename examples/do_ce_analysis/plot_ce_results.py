from nfoils.ce import CEPlotter

# select directory with calculated and experimental results files inside
# can analyse three different libraries: activities should be inside calc_results
experiment_directory = 'deuteron_nov24'
calc_results_file = 'c_results_fzk'
exp_results_file = 'e_results'

#input the plot labels for the three libraries you are analysing
# in same oder as in c_results file
libraries = ['TENDL-2021','IRDFF-II','ENDF/B-VIII']

# name for the C/E plot
plotname = 'test_fzk'

# set y axis size
y_upper = 2
y_lower = 0

# Flux normalisation for C results (if not already applied in fispact sim)
# calculated using be7 calcs from the xs_calculator program, but you do you
# 0.7(0.026985 rescale) or 0.801787(0.0309 rescale) for jendl/p-li. 0.8733 for initial unfolded175/p-li 
flux_norm_mean = 1
# Fractional uncertainty from the fux estimation
# 0.12783 for be7/p-li 
flux_frac_error = 0.12783 

#first weighted ave in list and last weighted ave in list for the WE calculations
first_we = 5
last_we = 20

# reorder the isotopes in your results files into capture-to-threshold
#order = [10,2,12,13, 7,5,6,9 ,4,3,0,1,14,11] #for proton_march_24
order = [10,2,15,16, 6,19,7,9,12,13,4,5,0,3,1,14,11] #for deuteron_nov_24

# how many isotopes (in the list above) are primarily thermal induced? 
# for split plot
plot_splitter = 6

#run 
ce_plotter = CEPlotter(experiment_directory,calc_results_file,exp_results_file,
                       plotname,y_upper,y_lower,flux_norm_mean,flux_frac_error,
                       first_we,last_we,libraries,order,plot_splitter)
ce_plotter.run()
