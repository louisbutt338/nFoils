""" 
module for doing flux estimation analysis on a lithium target

G Knoll, Radiation Detection and Measurement, 2010

M Majerle et al, Peak neutron production from the 7 Li(p,n) reaction 
in the 20-35 MeV range, EPJ Web of Conferences 239 (20010) 2020
"""

from math import pi, sqrt, log
import json


class TargetAnalysis:
    """ class to estimate flux from be7 activity in a lithium target
    """

    def __init__(self,target_json_name): 
        """ Initialise TargetAnalysis class

        Attributes
        ----------
        target_json_name : str
            name of the target data json
        """

        # set attributes
        self.target_json_name = target_json_name

        # load the target data
        with open(self.target_json_name) as target_file:
            self.target_data = json.load(target_file)

    def _no_of_isotopes(self,activity,t_half):
        """ calculate number of radioactive isotopes (knoll p2)

        Parameters
        ----------
        activity : float
            Activity of isotope
        t_half : float
            Halflife of isotope

        Returns
        -------
        isotopes : float
            No. of rad isotopes
        """
        activity_lambda = log(2)/t_half
        isotopes = activity/activity_lambda 
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

        Returns
        -------
        target_atoms : float
            No. of atoms in a target
        """
        avo_number = 6.02214076e23
        number_density = (avo_number*mass_density)/atom_mass
        target_atoms = number_density * thickness * pi * (radius)**2
        return target_atoms

    def _no_of_beam_particles(self,current,irrad_time):
        """ calculate number of charged particles incident on target
        for given current

        Parameters
        ----------
        current : float
            Current of charged particle in uA
        irrad_time : float
            Irradiation time in seconds

        Returns
        -------
        no_particles : float
            No. of particles in the beam
        """
        total_coulombs = current*1e-6*irrad_time
        no_particles = total_coulombs/(1.602176634e-19)
        return no_particles

    def _cross_section(self,no_isotopes,target_atoms,beam_flux):
        """ calculate a lithium 7 cross-section in millibarns (majerle)

        Parameters
        ----------
        no_isotopes : float
            Number of isotopes
        target_atoms: float
            Number of target atoms
        beam_flux : float
            Flux of the incident beam

        Returns
        -------
        lithium7_xs : float
            Real lithium 7 cross section 
        """
        lithium7_xs = (1e27*no_isotopes)/(target_atoms*beam_flux*0.925)
        return lithium7_xs

    def run(self,isotope_activity,isotope_halflife,current_list,
            timing_list,real_cross_section,current_rel_uncert,
            energy_rel_uncert):
        """ run the thing for be7 in a lithium target and print results

        Parameters
        ----------
        isotope_activity : list[float]
            Activity and relative uncertainty of the isotope in Bq
        isotope_halflife: float
            Halflife of the isotope
        current_list : list[float]
            List of faraday cup current readings in uA
        timing_list : list[int]
            List of irradiation timings in s
        real_cross_section : list[float]
            Cross section (mb) and relative uncertainty in a list
        current_frac_uncert : float
            Relative uncertainty on the current readings
        energy_frac_uncert : float
            Relative uncertainty on the energy of the incident
            particle

        Returns
        -------
        correction_factor : float
            correction factor to scale your Faraday cup current by
            to get current on the lithium target
        total_rel_uncert : float
            relative uncertainty on the above
        """

        # find target data from file 
        target_thickness = self.target_data["thickness"]
        target_mass_density = self.target_data["mass_density"]
        target_atomic_mass = self.target_data["atomic_mass"]
        target_radius = self.target_data["radius"]
        target_rel_uncert = self.target_data["relative_uncertainty"]
 
        # calculate incident particles and be7 cross section
        summed_incident_particles = sum(
            [self._no_of_beam_particles(
                j,timing_list[i]) for i,j in enumerate(current_list)])
        be7_cross_section = self._cross_section(
            self._no_of_isotopes(isotope_activity[0],
                                 isotope_halflife),
            self._no_of_target_atoms(target_thickness,
                                     target_mass_density,
                                     target_atomic_mass,
                                     target_radius),
            summed_incident_particles)
        
        # do calculations for the total relative uncertainty
        total_rel_uncert = sqrt(
            isotope_activity[1]**2
            +(real_cross_section[1])**2
            +target_rel_uncert**2
            +current_rel_uncert**2
            +energy_rel_uncert**2)

        # do calculation for the correction factor
        correction_factor = be7_cross_section/real_cross_section[0]
        print(f"flux correction factor from simulation is "
              f"{correction_factor} " 
              f"+- {100*total_rel_uncert} %")
        return correction_factor,total_rel_uncert