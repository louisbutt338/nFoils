""" example for plotting c/e results from experimental and calculated.
calculated results should have three datasets for three nuclear data libs,
 experimental should have one dataset from one experiment.
"""
from bfoils.ce import CEPlotter

#folder
folder = 'proton_march24'

# name for the C/E PNG plot
plotfile = 'bayes_unfold_fast_020326'

# initialise
ce_plotter = CEPlotter(folder,plotfile)

# name of calculated results and experimental results files
calc_results = 'c_results_bayes_unfold.json'
exp_results = 'e_results.json'

#input the plot labels for the libraries you are analysing
# in same oder as in c_results file
libraries = ['TENDL-2021','IRDFF-II','ENDF/B-VIII']

# Flux normalisation for C results (if not already applied in fispact sim)
# calculated using be7 calcs from the xs_calculator program
# 0.801787(used 0.0309 rescale) for model/p-li
# 1 (used proper rescale)for unfolded/pli 
# 1.1373 for model/dli nov24 results 
# (because you used the 0.1627225 rescale on the fispact calcs,
#  instead of 0.1847885)
# 1 for model/dli sep25
flux_norm_mean = 1
# Fractional uncertainty from the flux estimation
# 0.1285 for pli mar24
# 0.121225 for dli nov24
# 0.136018 for dli sep25
flux_frac_error =  0.1285

# reorder the isotopes in your results files into capture-to-threshold
#order = [10,2,12,13, 7,5,6,9 ,4,0,3,1,14,11]  # for pli mar24
order = [6,9 ,4,0,3,1,14,11]  # for just fast pli mar24
#order = [10,2,15,16, 6,19,7,9,12,4,5,0,3,  13,1,14,11]  # for dli nov24
#order = [8,2,13,14,  5,16,6,7,10,   4,0,3,11,1,12,9]  # for dli sep25

# first weighted ave in list and last weighted ave in new reordered list
# for the WE calculations e.g. [4,19] for isotopes 4-->19
we_isotopes = [0,19]

# which dataset would sir like to use for the weighted average calcs?
we_library = 'IRDFF-II'

# would sir like to add some vertical lines in between some isotopes? 
# i.e. for splitting the plot into thermal/threshold [4.5,13.5]
# .5 value will slot in between each isotope dataset:
# d1: [4.5,13.5] d2: [4.5,9.5] p:[5.5]
plot_split_list = [14.5]

# set y axis size e.g. [0.1,2.3] for 0-->2
y_axis = [0.1,3]

# set how far along the x axis to place the legend
legend_x_coord = 0.02 # d1/p1: 0.02, d2: 0.28

#run 
ce_plotter.run(calc_results,exp_results,libraries,flux_frac_error,
               we_isotopes,we_library,order,y_axis,flux_norm_mean,
               legend_x_coord,plot_split_list)

