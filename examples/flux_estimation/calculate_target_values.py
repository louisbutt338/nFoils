from nfoils.target import TargetAnalysis
import numpy as np

# be7 activity at end of irradiation and uncertainty (Bq), halflife of isotope
isotope_activity = [1e5,1e3]
isotope_halflife = 53.22*(24*3600)

# current array in uA, timing array in seconds
current_list = [8]   
timing_list  = [120] 

# target thickness (cm), radius (cm), material density (g/cm3), atomic mass (g/mol)
target_thickness = 0.1
target_radius = 0.1
target_mass_density = 0.1
target_atomic_mass = 5

# ENDFB8 7Li()7Be XS (mb) and fractional uncertainty 
# for the energy of particles in the middle of target
cross_section = [50,1]

# model neutron flux (n/cm2/src p) you want to rescale
foil_flux = 1e-05 

# model another neutron flux (n/cm2/src p) to rescale
target_flux = 0.005 

#run
analyse_target = TargetAnalysis(isotope_activity,isotope_halflife,current_list,
                                timing_list,target_thickness,target_radius,
                                target_mass_density,target_atomic_mass,
                                cross_section,foil_flux,target_flux)
analyse_target.run()