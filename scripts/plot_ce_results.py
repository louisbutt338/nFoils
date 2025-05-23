from nfoils.ce_plotter import CEPlotter

def main():
    # select working directory with 'calculated_activities' and 'experimental_activities' inside
    experiment = 'deuteron_nov24'
    working_directory = f'/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/analysis/{experiment}'
    calculated_results_file = 'calc_activities_60mb'
    experimental_results_file = 'exp_activities_apr25'

    # name for the C/E plot
    plotname = 'may2025_mcnp_TEST'

    # select analysis method out of root/interspec - should be included in the experimental_activities file
    experimental_analysis_method = 'interspec'

    # FLUX NORMALISATION for C results - calculated using be7 calcs from the xs_calculator program
    flux_norm_mean = 1 #0.7 for the mcnp/be7/p-li, 1 for unfolded175/be7/p-li 
    flux_percentage_error = 0.159 #0.159 for mcnp/be7/p-li , 0.218 for unfolded175/be7/p-li
    # redundant old normalisation factors
    #flux_norm_mean = 0.7 for the mcnp/be7/p-li, 1/0.4302 for unfolded175/be7/p-li 
    #flux_norm_error = 0.1496197
    #flux_norm_error = 0.5 * 0.040255/0.169795
    #flux_norm_mean = 0.6876

    #first weighted ave in list and last weighted ave in list for the WE calculations
    first_we = 0
    last_we = 5

    ce_plotter = CEPlotter(experiment,working_directory,calculated_results_file,
                           experimental_results_file,plotname,experimental_analysis_method,
                           flux_norm_mean,flux_percentage_error,first_we,last_we)
    ce_plotter.run()

if __name__ == '__main__':
    main() 