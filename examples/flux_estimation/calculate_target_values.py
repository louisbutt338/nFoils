""" example for doing flux estimation with lithium target
"""
from nfoils.target import TargetAnalysis

# path to json file with target data inside:
# target thickness (cm), radius (cm), material density (g/cm3),
# atomic mass (g/mol)
target_json = 'target_data'
analyse_target = TargetAnalysis(target_json)

# be7 activity at end of irradiation and uncertainty (Bq),
# halflife of isotope (s)
isotope_activity = [1e3,1e1]
isotope_halflife = 53.22*(24*3600)

# faraday cup current array in uA, timing array in seconds
current_list = [8]   
timing_list  = [120] 

# ENDFB8 7Li()7Be XS (mb) and fractional uncertainty 
# for the energy of particles in the middle of target
cross_section = [50,0.01]

# get correction factor for rescaling FC1 current values
# to target current values
c_factor,frac_uncert = analyse_target.run(isotope_activity,isotope_halflife,
                                          current_list,timing_list,
                                          cross_section)

# model neutron flux (n/cm2/src p) you want to rescale
foil_flux = 1e-05 

# calculate the rescaled flux
rescaled_flux = foil_flux*c_factor
print(f"rescaled flux is {rescaled_flux} n/cm2/src p +- {frac_uncert*100}%")
