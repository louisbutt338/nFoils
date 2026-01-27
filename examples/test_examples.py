""" test all examples: use to check package is installed correctly
"""

import subprocess

# test one example from input folder (str) and input script (str)
def test_example(folder,script):
    subprocess.run(['python',f'{script}.py'],
                    cwd=f'{folder}',
                    check=True,
                    stdout=subprocess.DEVNULL)

# test all examples here
print('testing examples')
test_example('analyse_fispact_json','get_fispact_activity')
test_example('analyse_gamma_spectrum','sum_ascii_files')
test_example('calibrate_gamma_detector','fit_efficiency_curve')
test_example('calibrate_gamma_detector','do_sum_correction')
test_example('estimate_target_flux','do_flux_estimation')
test_example('estimate_target_flux','do_rate_conversion')
test_example('extract_nuclear_data','get_response_functions')
test_example('extract_nuclear_data','get_spectrum_uncertainties')
test_example('calculate_foil_activities','do_activity_calculations')
test_example('calculate_foil_activities','do_ce_analysis')
test_example('calculate_foil_activities','get_unfolding_data')
test_example('unfold_neutron_spectrum','do_unfolding')
print('all nFoils examples run as expected - go measure some neutrons!')