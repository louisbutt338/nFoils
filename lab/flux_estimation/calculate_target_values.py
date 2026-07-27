""" example for doing flux estimation with lithium target. 
calculates conversion for FC1 to target, charged particle rate,
and writes a fispact history
"""
import json
import numpy as np
from bfoils.target import TargetAnalysis
from bfoils.history import IrradTimeline

# calculating the correction factor from cup current to target current

# path to json file with target data inside:
# target thickness (cm), radius (cm), material density (g/cm3),
# atomic mass (g/mol) and their total relative uncertainty
target_json = 'target_data.json'
analyse_target = TargetAnalysis(target_json)

# activity at end of irradiation and relative uncertainty (Bq)
# p:[1.5979e5,0.04448] (new cal) 
# d1: [2.15488e5,0.04469] or [2.15170e+05,0.04462]
# d2: [5.44243e+05, 0.04445]
isotope_activity = [2.15170e+05,0.04462]
# halflife of be7 isotope
isotope_halflife = 53.22*(24*3600)

# current array in uA, timing array in seconds
# grab deuteron experiment uA data from a json or just list here
with open('deuteron_sep25/cup_currents.json'  # d1:'deuteron_nov24...'
          ) as cup_json:
    cup_data = json.load(cup_json)
current_list = cup_data["cup_currents_ua"]
timing_list = cup_data["cup_timings_s"]
# p:
#current_list = [5.5,9,10]    # p:[5.5,9,10] or [8.778]
#timing_list  = [i*60 for i in [20,67,41.5]]  # p:[20,67,41.5] or [128.5]

# relative uncertainty on the current (0.06)
current_rel_uncert = 0

# ENDFB8 7Li()7Be XS (mb) and fractional uncertainty
# for the energy of particles in the middle of target
# p:[33.5,0.0324] d1:[64.3,0.05] d2:[50,0.05]
cross_section = [64.3,0.05]

# relative uncertainty on the energy of the incident particle
# p: 0.0125 d1: 0.0143 d2: 0.02062
energy_rel_uncert = 0.0143

# get correction factor for rescaling FC1 current values 
# to target current values
correction,rel_uncert = analyse_target.run(isotope_activity,isotope_halflife,
                                           current_list,timing_list,
                                           cross_section,current_rel_uncert,
                                           energy_rel_uncert)


# estimating target source strength and foil neutron flux

# mcnp model fe foil neutron flux total (n/cm2/src p)
foil_flux = 5.49226e-5  # p: 1.36564e-5 d1: 5.49226e-5 d2: 3.45995e-5
# mcnp model li target neutron flux total (n/cm2/src p)
target_flux = 4.9867e-3 # p:2.7884e-3 d1: 4.9867e-3 d2: 2.68287e-3

# get target area for source strength:
with open(f'{target_json}.json'
          ) as target_file:
    target_data = json.load(target_file)
target_area = np.pi*(target_data["radius"]**2)

# input MCNP value for number of charged particles per s per uA
source_p_per_s_per_ua = 6.24151e+12

# calculate source strength and flux
source_strength = (correction*source_p_per_s_per_ua*target_area*target_flux)
flux = correction*source_p_per_s_per_ua*foil_flux
print("at FC1 = 1uA:")
print(f"source strength is {source_strength:.3e} "
      f"+- {rel_uncert*source_strength:.3e} n/s")
print(f"flux is {flux:.3e} "
      f"+- {rel_uncert*flux:.3e} n/cm2/s")


# write a fispact history for particles/s on the target

# calculate target currents so fispact calculations are correct
target_current_list = [i*correction for i in current_list]

# set fispact irrad hist filename
output_filename = "deuteron_nov24/fispact_history_target"

#initialise and run
history = IrradTimeline(target_current_list,timing_list)
history.fispact_hist_writer(source_p_per_s_per_ua,output_filename)
