from nfoils.target import TargetAnalysis
import json
import numpy as np

#grab deuteron experiment uA data
deuteron_currents_data = json.load(open('deuteron_experiment_currents.json'))

# activity at end of irradiation and uncertainty (Bq), halflife of isotope
isotope_activity = [189476,44921] #p:[1.5979e5,7.10764e3] (new cal) d:[189476,44921] 
isotope_halflife = 53.22*(24*3600)

# DEUTERONS: current array in uA, timing array in seconds
current_list = deuteron_currents_data["currents_in_ua"]
timing_list = np.ones(5530)
# PROTONS: current array in uA, timing array in seconds
#current_list = [5.5,9,10]    # p:[5.5,9,10] or [8.778]
#timing_list  = [i*60 for i in [20,67,41.5]]  # p:[20,67,41.5] or [128.5]

# target thickness (cm), radius (cm), material density (g/cm3), atomic mass (g/mol)
target_thickness = 0.05
target_radius = 0.60
target_mass_density = 0.534
target_atomic_mass = 6.941

# ENDFB8 7Li()7Be XS (mb) and fractional uncertainty for the energy of particles in the middle of target
cross_section = [64.3,0.0324] # p:33.5 d:64.3

# model neutron flux (n/cm2/src p) to rescale
foil_flux = 5.49226e-05 # p:1.36564e-05 d:5.49226e-05

# model li target neutron flux (n/cm2/src p) to rescale src strength
target_flux = 0.004986699 # p:2.788406e-3 d: 0.004986699

#run
analyse_target = TargetAnalysis(isotope_activity,isotope_halflife,current_list,
                                timing_list,target_thickness,target_radius,
                                target_mass_density,target_atomic_mass,
                                cross_section,foil_flux,target_flux)
analyse_target.run()