""" example for doing flux estimation with lithium target.
Compares your current readings to expected cross section value
and outputs a correction factor to give you current on the target
"""
from nfoils.target import TargetAnalysis

# path to json file with target data inside:
# target thickness (cm), radius (cm), material density (g/cm3),
# atomic mass (g/mol) and their total relative uncertainty
target_json = 'target_data'
analyse_target = TargetAnalysis(target_json)

# be7 activity at end of irradiation and relative uncertainty (Bq),
# halflife of isotope (s)
isotope_activity = [1e3,0.01]
isotope_halflife = 53.22*(24*3600)

# faraday cup current array in uA, timing array in seconds
# and relative uncertainty on the current
current_list = [8]
timing_list  = [120]
current_rel_uncert = 0.06

# ENDFB8 7Li()7Be XS (mb) and relative uncertainty 
# for the energy of particles in the middle of target
cross_section = [50,0.01]

# relative uncertainty on the enegy of the incident particle
energy_rel_uncert = 0.0125

# get correction factor for rescaling FC1 current values
# to target current values
c_factor,rel_uncert = analyse_target.run(isotope_activity,isotope_halflife,
                                         current_list,timing_list,
                                         cross_section,current_rel_uncert,
                                         energy_rel_uncert)

# model neutron flux (n/cm2/src p) you want to rescale
foil_flux = 1e-05 

# calculate the rescaled flux
rescaled_flux = foil_flux*c_factor
print(f"rescaled flux is {rescaled_flux} n/cm2/src p +- {rel_uncert*100}%")
