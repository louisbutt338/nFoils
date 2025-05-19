from math import pi, sqrt, log, exp
import numpy as np # type: ignore

class TargetAnalysis:
    def __init__(self, isotope_activity,isotope_halflife,current_list,timing_list,
               target_thickness,target_radius,target_mass_density,target_atomic_mass,
               real_cross_section,real_source_strength,real_fe_flux):
        self.isotope_activity = isotope_activity
        self.isotope_halflife = isotope_halflife
        self.current_list =  current_list
        self.timing_list = timing_list
        self.target_thickness = target_thickness
        self.target_radius = target_radius
        self.target_mass_density = target_mass_density
        self.target_atomic_mass = target_atomic_mass
        self.real_cross_section = real_cross_section
        self.real_source_strength = real_source_strength
        self.real_fe_flux = real_fe_flux

    # calculate number of be7 atoms
    def _no_of_isotopes(self,activity,t_half):
        activity_lambda = log(2)/t_half
        return activity/activity_lambda #activity/(1-exp(-irrad_time*activity_lambda))

    # calculate number of target atoms
    def _no_of_target_atoms(self,thickness,mass_density,atom_mass,radius):
        avo_number = 6.02214076e23
        number_density = (avo_number*mass_density)/atom_mass
        return number_density * thickness * np.pi * (radius)**2

    # calculate charged particles incident on target for each current
    def _no_of_beam_particles(self,current,irrad_time):
        total_coulombs = current*1e-6*irrad_time
        no_particles = total_coulombs/(1.602176634e-19)

        # if particle == 'proton':
        #     scaling_factor_1ua = 6.24151e12
        # if particle == 'deuteron':
        #     scaling_factor_1ua = 6.24151e13
        # particle_flux = total_flux*scaling_factor_1ua*current

        return no_particles 

    # calculate cross-section in millibarns
    def _cross_section(self,no_isotopes,target_atoms,beam_flux):
        return (1e27*no_isotopes)/(target_atoms*beam_flux*0.925)

    def run(self):
        # do calculations for the fractional uncertainty 
        summed_incident_particles = sum([self._no_of_beam_particles(self.current_list[i],self.timing_list[i]) for i in range(len(self.current_list))])
        be7_cross_section = self._cross_section(self._no_of_isotopes(self.isotope_activity[0],self.isotope_halflife),self._no_of_target_atoms(self.target_thickness,self.target_mass_density,self.target_atomic_mass,self.target_radius),summed_incident_particles)
        total_frac_uncert = np.sqrt((self.isotope_activity[1]/self.isotope_activity[0])**2 + (self.real_cross_section[1]/self.real_cross_section[0])**2)
        print(f"total flux-estimation fractional uncert = {total_frac_uncert} " )

        # do calculation for the factor the XS is wrong by (for the timing/current arrays used here) - use to estimate #neutrons per second from lithium and neutron flux on iron foil
        correction_factor = be7_cross_section/self.real_cross_section[0]
        print(f"flux correction factor from simulation is {correction_factor} +- {correction_factor*total_frac_uncert}")

        # do calculations for experimental source strength of lithium target - by benchmarking to mcnp value
        source_particles_per_second_10ua = 6.24151e+13
        source_strength = correction_factor*source_particles_per_second_10ua*self.real_source_strength
        flux = correction_factor*source_particles_per_second_10ua*self.real_fe_flux
        print(f"Source strength at 10uA is {source_strength:.5e} +- {total_frac_uncert*source_strength:.5e} n/s")
        print(f"Flux at fe foil at 10uA is {flux:.5e} +- {total_frac_uncert*flux:.5e} n/cm2/s")


