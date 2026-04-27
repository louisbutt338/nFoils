""" example for doing flux estimation with lithium target. 
calculates conversion for FC1 to target, charged particle rate,
and writes a fispact history
"""
import json
import numpy as np
from nfoils.target import TargetAnalysis
from nfoils.history import IrradTimeline

# calculating the correction factor from cup current to target current

# path to json file with target data inside:
# target thickness (cm), radius (cm), material density (g/cm3),
# atomic mass (g/mol) and their total relative uncertainty
target_json = 'target_data.json'
analyse_target = TargetAnalysis(target_json)

# activity at end of irradiation and relative uncertainty (Bq)
# p:[1.5979e5,0.04448] (new cal) 
isotope_activity = [1.5979e5,0.04448]
# halflife of be7 isotope
isotope_halflife = 53.22*(24*3600)

# current array in uA, timing array in seconds
current_list = [5.5,9,10]
timing_list  = [i*60 for i in [20,67,41.5]] 

# relative uncertainty on the current (0.06)
current_rel_uncert = 0

# ENDFB8 7Li()7Be XS (mb) and fractional uncertainty
# for the energy of particles in the middle of target
cross_section = [33.5,0.0324]

# relative uncertainty on the energy of the incident particle
energy_rel_uncert = 0.0125

# get correction factor for rescaling FC1 current values 
# to target current values
correction,rel_uncert = analyse_target.run(isotope_activity,isotope_halflife,
                                           current_list,timing_list,
                                           cross_section,current_rel_uncert,
                                           energy_rel_uncert)


# estimating target source strength and foil neutron flux

# mcnp model fe foil neutron flux total (n/cm2/src p)
foil_flux = 1.36564e-5  
# mcnp model li target neutron flux total (n/cm2/src p)
target_flux = 2.7884e-3 

# get target area for source strength:
with open(f'{target_json}'
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
output_filename = "proton_march24/fispact_history_target"

#initialise and run
history = IrradTimeline(target_current_list,timing_list)
history.fispact_hist_writer(source_p_per_s_per_ua,output_filename)
