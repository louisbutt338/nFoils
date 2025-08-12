from nfoils.ce import CEPlotter

# select directory with calculated and experimental results files inside
experiment_directory = '../do_ce_analysis'

# name for the C/E plot
plotname = 'test_fzk'

# set y axis size
y_upper = 2
y_lower = 0

# name of calc and exp results files
calc_results_file = 'c_results'
exp_results_file = 'e_results'

#input the plot labels for the three libraries you are analysing
# in same oder as in c_results file
libraries = ['TENDL-2021','IRDFF-II','ENDF/B-VIII']

# Flux normalisation to scale C results
flux_norm_mean = 1.1

# Fractional uncertainty on the fux estimation
flux_frac_error = 0.05

#first weighted ave in list and last weighted ave in list for the WE calculations
first_we = 5
last_we = 20

# reorder the isotopes in your results files into capture-to-threshold
order = [10,2,15,16, 6,19,7,9,12,13,4,5,0,3,1,14,11] 

# how many isotopes (in the list above) are primarily thermal induced? 
# for splitting plot
plot_splitter = 4

#run 
ce_plotter = CEPlotter(experiment_directory,calc_results_file,exp_results_file,
                       plotname,y_upper,y_lower,flux_norm_mean,flux_frac_error,
                       first_we,last_we,libraries,order,plot_splitter)
ce_plotter.run()
