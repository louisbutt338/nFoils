""" example for plotting c/e results from experimental and calculated

run after 1_calculate_exp_activities.py
"""
from nfoils.ce import CEPlotter

# path to the calculated and experimental results files
folder = 'example'

# name to save the C/E plot under as PNG (saves in same folder)
plotname = 'test'

# initialise object
ce_plotter = CEPlotter(folder,plotname)

# names of calculated and experimental json results files
calc_results = 'c_results'
exp_results = 'e_results'

#input the plot labels for the three libraries you are analysing
# in same oder as in c_results file
libraries = ['TENDL-2021','IRDFF-II','ENDF/B-VIII']

# Flux normalisation to scale C results
flux_norm_mean = 1.1

# Fractional uncertainty on the fux estimation
flux_frac_error = 0.05

# reorder the isotopes in your results files into capture-to-threshold
order = [10,2,15,16, 6,19,7,9,12,13,4,5,0,3,1,14,11] 

# first weighted ave in list and last weighted ave in new reordered list
# for the WE calculations e.g. [4,19] for isotopes 4-->19
we_isotopes = [5,20]

# which dataset would sir like to use for the weighted average calcs?
we_library = 'ENDF/B-VIII'

# how many isotopes (in the list above) are primarily thermal induced? 
# for splitting plot
plot_splitter = 4

# set y axis size e.g. [0,2] for 0-->2
y_axis = [0.1,2.3]

#do c/e plot and weighted ave analysis
ce_plotter.run(calc_results,exp_results,flux_norm_mean,
               flux_frac_error,we_isotopes,libraries,
               order,plot_splitter,y_axis, we_library)
