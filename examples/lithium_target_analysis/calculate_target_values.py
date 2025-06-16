from nfoils.xs_calculator import TargetAnalysis
import json

#grab deuteron experiment uA data
deuteron_currents_file = 'deuteron_experiment_currents.json'
deuteron_currents_data = json.load(open(deuteron_currents_file))

# activity at end of irradiation and uncertainty (Bq), halflife of isotope
isotope_activity = [1.5979e5,7.10764e3]    #p:[1.5979e5,7.10764e3] (new cal) d:[189476,44921] 
isotope_halflife = 53.22*(24*3600)

# DEUTERONS: current array in uA, timing array in seconds
#current_list = deuteron_currents_data["currents_in_ua"]
#timing_list = np.ones(5530)
# PROTONS: current array in uA, timing array in seconds
current_list = [5.5,9,10]    # p:[5.5,9,10]   or [8.778]  d:[10]
timing_list  = [i*60 for i in [20,67,41.5]]  # p:[20,67,41.5] or [128.5]  d:[20]

# target thickness (cm), radius (cm), material density (g/cm3), atomic mass (g/mol)
target_thickness = 0.05
target_radius = 0.60
target_mass_density = 0.534
target_atomic_mass = 6.941

# real 7Li XS (mb) for the ave energy of particles in the target to get correction factor from
# and real target source strength (n per src p) and flux (per source proton) on Fe foil 
real_cross_section = [33.5,1.0857] # p:33 Using hand-drawn to match with ENDFB8/experimental data  #d:60 or 85 ? Two hand-drawn lines as two datasets
real_source_strength = 3.28563E-03 # p:3.28563E-03 d: 
real_fe_flux = 1.37267e-05 # p:1.37267e-05 d:3.42800e+09/6.24151e+13

#run
analyse_target = TargetAnalysis(isotope_activity,isotope_halflife,current_list,timing_list,
                                    target_thickness,target_radius,target_mass_density,target_atomic_mass,
                                    real_cross_section,real_source_strength,real_fe_flux)
analyse_target.run()