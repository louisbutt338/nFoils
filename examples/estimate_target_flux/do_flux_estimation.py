""" do a flux estimation for lithium target neutron source. 
(1) Compare your faraday cup current readings to expected cross section value
and output a correction factor to give you real current on the target 
(2) Rescale an input value by this correction i.e. neutron flux 
(3) Print a fispact irradiation history for the target
"""
from nfoils.target import TargetAnalysis
from nfoils.history import IrradTimeline

# calculate the correction factor and uncertainty:

# set path to json file with target data inside and initialise
# should include single dictionary including:
#  target thickness (cm)
#  radius (cm)
#  material density (g/cm3)
#  atomic mass (g/mol)
#  total relative uncertainty on all the above
target_json = 'data/target'
analyse_target = TargetAnalysis(target_json)

# be7 activity (Bq) at end of irradiation and relative uncertainty
isotope_activity = [1e3,0.01]
# halflife of be7 (s)
isotope_halflife = 53.22*(24*3600)

# irradiation experiment timeline details
# faraday cup current array in uA
current_list = [8]
# relative uncertainty on the current
current_rel_uncert = 0.06
# timing array in seconds
timing_list  = [120]

# cross-section for 7Li()7Be in mb and relative uncertainty 
# for the particle energy in the middle of target
cross_section = [50,0.01]

# relative uncertainty on the energy of the incident particle
energy_rel_uncert = 0.0125

# get correction factor to rescale cup current values to target current values
c_factor,rel_uncert = analyse_target.run(isotope_activity,isotope_halflife,
                                         current_list,timing_list,
                                         cross_section,current_rel_uncert,
                                         energy_rel_uncert)


# example for rescaling an input value by the correction:

# example value you want to rescale i.e. neutron flux (n/cm2/src p)
foil_flux = 1e-05 

# input MCNP value for number of charged particles per s per uA
source_p_per_s_per_ua = 6.24151e+12

# calculate the rescaled neutron flux
rescaled_flux = foil_flux*c_factor*source_p_per_s_per_ua
print(f"rescaled flux is {rescaled_flux} n/cm2/s/uA +- {rel_uncert*100}%")


# write a fispact history in units of particles/s on the target:

# rescale cup currents to target currents
target_current_list = [i*c_factor for i in current_list]

# set fispact irrad hist filename
output_filename = "fispact_history_target"

# convert target currents into target rate (particles/s) and print 
# a fispact history with them in 
history = IrradTimeline(target_current_list,timing_list)
history.fispact_hist_writer(source_p_per_s_per_ua,output_filename)
