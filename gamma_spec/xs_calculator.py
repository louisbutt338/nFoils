from math import pi, sqrt, log, exp
import numpy as np # type: ignore

################################################################################
########### USER INPUTS ########################################################
################################################################################
################################################################################

# general XS calculator for proton/deuteron beam on a THIN target

# activity at end of irradiation and uncertainty (Bq)
isotope_activity = [189476,44921] #[202151,43070]  #[189476,44921]
isotope_halflife = 53.22*(24*3600)

# array of currents in uA, array of timings in minutes
current_list = [10] #[5.5,9,10]   #[10]#[8.778]#
timing_list  = [20] #[20,67,41.5] #[20]#[128.5]#

# target thickness (cm) and target radius (cm)
target_thickness = 0.05
target_radius = 0.6

# target material density (g/cm3) and atomic mass (g/mol)
target_mass_density = 0.534
target_atomic_mass = 6.941

###########################################################################################
###########################################################################################
###########################################################################################
##########################################################################################

timing_list_seconds = [i*60 for i in timing_list]
print(timing_list_seconds)

# calculate number of be7 atoms
def no_of_isotopes(activity,t_half):
    activity_lambda = log(2)/t_half
    return activity/activity_lambda #activity/(1-exp(-irrad_time*activity_lambda))

# calculate number of target atoms
def no_of_target_atoms(thickness,mass_density,atom_mass,radius):
    avo_number = 6.02214076e23
    number_density = (avo_number*mass_density)/atom_mass
    return number_density * thickness * np.pi * (radius)**2

# calculate charged particles incident on target for each current
def no_of_beam_particles(current,irrad_time):
    total_coulombs = current*1e-6*irrad_time
    no_particles = total_coulombs/(1.602176634e-19)

    # if particle == 'proton':
    #     scaling_factor_1ua = 6.24151e12
    # if particle == 'deuteron':
    #     scaling_factor_1ua = 6.24151e13
    # particle_flux = total_flux*scaling_factor_1ua*current

    return no_particles 

# calculate cross-section in millibarns
def cross_section(no_isotopes,target_atoms,beam_flux):
    return (1e27*no_isotopes)/(target_atoms*beam_flux*0.95)

print(f"Be7 XS = {cross_section(no_of_isotopes(isotope_activity[0],isotope_halflife),no_of_target_atoms(target_thickness,target_mass_density,target_atomic_mass,target_radius),no_of_beam_particles(current_list[0],timing_list_seconds[0]))} mb" )

# adds together from diff irradiation periods

