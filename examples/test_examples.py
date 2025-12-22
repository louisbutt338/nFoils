""" test functionality of all examples. 
use to check package is working as intended
"""

import subprocess

# function to test each example
def test_example(folder,script):
    """ test a single example from input folder (str) and input script (str)
    """
    subprocess.run(['python',f'{script}.py'],
                    cwd=f'{folder}',
                    check=True)

# test all examples here
test_example('fispact_jsons','postprocess_json')
test_example('nuclear_data','get_library_data')
test_example('nuclear_data','get_spectrum_uncertainties')
test_example('activity_analysis','calculate_e_results')
test_example('activity_analysis','do_ce_analysis')
test_example('activity_analysis','get_reaction_rates')
test_example('flux_estimation','calculate_target_values')
test_example('flux_estimation','rate_conversion')
test_example('gamma_spec','fit_efficiency_curve')
test_example('gamma_spec','sum_correct_source')
test_example('ascii_files','sum_and_plot_asciis')
test_example('unfolding','unfold_spectrum')
print('***','Congrats! All nFoils examples run as expected','***')