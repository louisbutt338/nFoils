""" 
module for doing analysis on a lithium target
"""

from math import pi, sqrt, log
import numpy as np 

class TargetAnalysis:
    """ class to estimate flux from be7 activity
    """
    def __init__(self, isotope_activity,isotope_halflife,
                 current_list,timing_list,target_thickness,
                 target_radius,target_mass_density,
                 target_atomic_mass, real_cross_section,
                 real_foil_flux,real_target_flux):
        """ Initialise class
        """

        # set the parameters
        self.isotope_activity = isotope_activity
        self.isotope_halflife = isotope_halflife
        self.current_list =  current_list
        self.timing_list = timing_list
        self.target_thickness = target_thickness
        self.target_radius = target_radius
        self.target_mass_density = target_mass_density
        self.target_atomic_mass = target_atomic_mass
        self.real_cross_section = real_cross_section
        self.real_foil_flux = real_foil_flux
        self.real_target_flux = real_target_flux

    def _no_of_isotopes(self,activity,t_half):
        """ calculate number of atoms

        Parameters
        ----------
        activity : float
            Activity of isotope
        t_half : float
            Halflife of isotope
        """
        activity_lambda = log(2)/t_half
        isotopes = activity/activity_lambda 
        #isotopes = activity/(1-exp(-irrad_time*activity_lambda))
        return isotopes 

    def _no_of_target_atoms(self,thickness,mass_density,atom_mass,radius):
        """ calculate number of atoms in a target

        Parameters
        ----------
        thickness : float
            Target thickness in cm
        mass_density : float
            Density of target material in g/cm3
        atom_mass : float
            Atomic mass of the target material
        radius : float
            Radius of the target in cm
        """
        avo_number = 6.02214076e23
        number_density = (avo_number*mass_density)/atom_mass
        return number_density * thickness * pi * (radius)**2

    def _no_of_beam_particles(self,current,irrad_time):
        """ calculate number of charged particles
          incident on target for given current

        Parameters
        ----------
        current : float
            Current of charged particle in uA?
        irrad_time : float
            Irradiation time in seconds
        """
        total_coulombs = current*1e-6*irrad_time
        no_particles = total_coulombs/(1.602176634e-19)
        # if particle == 'proton':
        #     scaling_factor_1ua = 6.24151e12
        # if particle == 'deuteron':
        #     scaling_factor_1ua = 6.24151e13
        # particle_flux = total_flux*scaling_factor_1ua*current

        return no_particles 

    def _cross_section(self,no_isotopes,target_atoms,beam_flux):
        """ calculate a cross-section in millibarns

        Parameters
        ----------
        no_isotopes : float
            Number of isotopes
        target_atoms: float
            Number of target atoms
        beam_flux : float
            Flux of the incident beam
        """
        return (1e27*no_isotopes)/(target_atoms*beam_flux*0.925)

    def run(self):
        """ run the thing for be7 in a lithium target
        """
        # do calculations for the fractional uncertainty 
        summed_incident_particles = sum(
            [self._no_of_beam_particles(
                self.current_list[i],
                self.timing_list[i]) for i in range(
                    len(self.current_list))])
        be7_cross_section = self._cross_section(
            self._no_of_isotopes(self.isotope_activity[0],
                                 self.isotope_halflife),
            self._no_of_target_atoms(self.target_thickness,
                                     self.target_mass_density,
                                      self.target_atomic_mass,
                                      self.target_radius),
                                      summed_incident_particles)
        total_frac_uncert = sqrt(
            (self.isotope_activity[1]/self.isotope_activity[0])**2
            +(self.real_cross_section[1])**2)
        print(f"flux-estimation fractional uncert (without FC1 uncert and incident proton energy uncert)= {total_frac_uncert} " )

        # do calculation for the correction factor
        correction_factor = be7_cross_section/self.real_cross_section[0]
        print(f"flux correction factor from simulation is {correction_factor} " 
              f"+- {correction_factor*total_frac_uncert}")

        # do estimations for source strength of target and iron foil neutron flux
        # by benchmarking to mcnp values
        source_p_per_s_10ua = 6.24151e+13
        target_area = np.pi*(self.target_radius**2)
        source_strength = (correction_factor*source_p_per_s_10ua*
                           target_area*self.real_target_flux)
        flux = correction_factor*source_p_per_s_10ua*self.real_foil_flux
        print(f"rescaled source strength when FC1=10uA is {source_strength:.5e} "
              f"+- {total_frac_uncert*source_strength:.5e} n/s")
        print(f"rescaled flux when FC1=10uA is {flux:.5e} "
              f"+- {total_frac_uncert*flux:.5e} n/cm2/s")


