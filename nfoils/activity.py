"""
module for calculating activity for a whole ton of isotopes

G Knoll, Radiation Detection and Measurement, 2010
D Arnold et al, Fundamentals of gamma spectrometry, 2018
"""
import json
from datetime import datetime
from math import pi, sqrt, log, exp
import numpy as np
import actigamma as ag
from scipy.integrate import quad


class ActivityCalc:
    """ class that calculates activities for isotope data written in a json
    """
    def __init__(self, data_file_name, json_path, irradiation_end):
        """ Initialise class
        """

        # set the parameters
        self.data_file_name = data_file_name
        self.json_path = json_path
        self.irradiation_end = irradiation_end

        # load the data
        with open(f'{self.json_path}/{self.data_file_name}.json'
                  ) as json_datafile:
            self.json_file_data = json.load(json_datafile)

    def _decay_time(self, isotope_name):
        """ calculating the decay time
        from input irradiation end and measurement start

        Parameters
        ----------
        isotope_name : str
            Name of isotope being run

        Returns
        -------
        dec_time : float
            Decay time of the isotope
        """
        dt_elements = self.json_file_data[isotope_name]['datetime']
        ts = datetime(*dt_elements)
        dec_time = (ts-self.irradiation_end).total_seconds()
        return dec_time

    def _get_decay_database(self, isotope_name):
        """ read the pypact actigamma decay2012 database
        and get halflife, intensities,peaks for specified isotope

        Parameters
        ----------
        isotope_name : str
            Name of isotope being run

        Returns
        -------
        i : list[float]
            Array of gamma peak intensities
        e : list[float]
            Array of associated gamma energies
        half_life : float
            Half life of the input isotope
        """
        SPECTYPE = "gamma"
        db = ag.Decay2012Database()
        half_life = db.gethalflife(isotope_name)
        intensity = db.getintensities(isotope_name, spectype=SPECTYPE)
        peak_energy_kev = 1e-3 * db.getenergies(isotope_name,
                                                spectype=SPECTYPE)
        sorted_lists = sorted(zip(peak_energy_kev, intensity),
                              key=lambda x: x[1])
        e, i = zip(*reversed(sorted_lists))
        return (i, e, half_life)

    def _solid_angle(self, crystal_radius, distance):
        """ solid angle approximation (knoll, p120)

        Parameters
        ----------
        crystal_radius : float
            Radius of detector crystal in cm
        distance : float
            Distance from detector to foil in cm

        Returns
        -------
        solid_ang : float
            the solid angle for point source
        """
        solid_ang = 2 * pi * (1 - distance / sqrt(distance**2
                                                  + crystal_radius**2))
        return solid_ang

    def _solid_angle_disc(self, crystal_radius, distance, foil_radius):
        """ solid angle approximation for foils/discs (knoll, p121)

        Parameters
        ----------
        crystal_radius : float
            Radius of detector crystal in cm
        distance : float
            Distance from detector to foil in cm
        foil_radius : float
            Radius of the foil in cm

        Returns
        -------
        solid_ang_disc : float
            the solid angle for a disc source
        """
        alpha = (foil_radius / distance)**2
        beta = (crystal_radius / distance)**2
        f_1 = (((5/16) * (beta/((1+beta)**(7/2))))
               - ((35/64) * ((beta**2)/((1+beta)**(9/2)))))
        f_2 = (((35/128) * (beta/((1+beta)**(9/2))))
               - ((315/256) * ((beta**2)/((1+beta)**(11/2))))
               + ((1155/1024) * ((beta**3)/((1+beta)**(13/2)))))
        sol_ang_disc = (2 * pi * (1 - (1/((1+beta)**(1/2)))
                        - (3/8)*((alpha*beta)/((1+beta)**(5/2)))
                        + (alpha**2)*f_1 - (alpha**3)*f_2))
        return sol_ang_disc

    def _efficiency_abs(self, energy, n0, n1, n2, n3):
        """ equation for the log-polynomial efficiency curves (knoll, p458)

        Parameters
        ----------
        energy : float
            Gamma energy in keV
        n0 : float
            Zero term parameter
        n1 : float
            1st term parameter for E^1
        n2 : float
            1st term parameter for E^2
        n3 : float
            1st term parameter for E^3

        Returns
        -------
        eff : float
            the efficiency of the detector
        """
        polynomial = (n0 + n1*np.log(energy)**1
                      + n2*np.log(energy)**2
                      + n3*np.log(energy)**3)
        eff = np.exp(polynomial)
        return eff

    def _activity_livetime(self, c, i, e, foil_dist, isotope_name):
        """ use the foil measurement distance and efficiency curves
        to calculate activity over the live time from counts 
        (knoll, p120) (arnold, p42)

        Parameters
        ----------
        c : float
            Counts
        i : float
            Intensity of gamma peak
        e : float
            Energy in keV of gamma peak
        foil_dist : float
            Distance from detector to foil in cm
        isotope_name : str
            Name of isotope

        Returns
        -------
        activity : float
            the foil activity over detector livetime
        """
        # b03 hpge values
        if foil_dist == 1:
            detector_radius = 3.25
            eff_values = [-23.385, 11.3457, -1.9302, 0.10197]
        # g11 bege values
        if foil_dist == 0.4:
            detector_radius = 3.75
            eff_values = [-25.7209, 13.8415, -2.48, 0.1352]
        # b03 hpge values
        if foil_dist == 38:
            detector_radius = 3.25
            eff_values = [-27.773, 11.13964, -1.8975, 0.10025]
        foil_rad = self.json_file_data[isotope_name]['foil_radius_cm']
        solid_angle_ratio = (self._solid_angle_disc(detector_radius,
                                                    foil_dist, foil_rad)
                             / self._solid_angle(detector_radius, foil_dist))
        selected_efficiency = (self._efficiency_abs(e, *eff_values)
                               * solid_angle_ratio)
        activity = c / (i * selected_efficiency)
        return activity

    def _self_attenuation_correction(self, material, e, thickness, density):
        """ gamma self absorption correction factor taken from XCOM mu data
        (arnold, p41)

        Parameters
        ----------
        material : str
            Foil material
        E : float
            Energy (keV) of gamma peak
        thickness : float
            Thickness of foil in cm
        density : float
            Density of foil material in g/cm3

        Returns
        -------
        self_att_factor : float
            self attenuation factor for the foil
        """
        xcom = np.fromfile(f'../../data/XCOM_new/{material}.txt', sep=" ")
        mass_coeff = np.interp(e/1000, xcom[::2], xcom[1::2])
        self_att_factor = ((mass_coeff * density * thickness)
                           / (1 - exp(- mass_coeff * density * thickness)))
        return self_att_factor

    def _activity_integrand(self, t, half_life):
        """ integral shortcut for activity interpolation (arnold, p44)

        Parameters
        ----------
        t : float
            time (s)
        half_life : float
            Half-life of isotope (s)

        Returns
        -------
        intergrand : float
            activity integral from time and halflife
        """
        integrand = exp(- log(2) * (t/half_life))
        return integrand

    def _activity_0(self, c, i, e, measurement_distance,
                    isotope_name, halflife):
        """ calculate initial activity after irrad (w/o corrections)
        from the measured activity over a live time (arnold, p44)

        Parameters
        ----------
        c : float
            Counts
        i : float
            Intensity of gamma peak
        e : float
            Energy in keV of gamma peak
        measurement_distance : float
            Distance from detector to foil in cm
        isotope_name : str
            Name of isotope
        halflife : float
            halflife of isotope

        Returns
        -------
        activity : float
            uncorrected activity at the end of irradiation
        """
        decay_time = self._decay_time(isotope_name)
        top_time_band = (decay_time
                         + self.json_file_data[isotope_name]['live_time'])
        quad_integral = quad(self._activity_integrand,
                             decay_time, top_time_band,
                             args=halflife)
        activity = (self._activity_livetime(c, i, e, measurement_distance,
                                            isotope_name)
                    / (quad_integral))[0]
        return activity

    def _reaction_rates(self, a, irradiation_time, halflife):
        """ calculate ave reaction rate
        from the activity a0 at time t0 under irradiation

        Parameters
        ----------
        a : float
            Activity after irradiaiton
        irradiation_time : float
            Time the irradiation spanned in s
        halflife : float
            Halflife of isotope in s

        Returns
        -------
        rr_ave : float
            average reaction rate for the irradiation
        """
        rr_ave = a / (1 - self._activity_integrand(irradiation_time, halflife))
        return rr_ave

    def _run_one_isotope(self, isotope_name, measurement_distance,
                         eff_uncert, irrad_time):
        """ run analysis for one isotope

        Parameters
        ----------
        measurement_distance : float
            Distance from detector to foil in cm
        isotope_name : str
            Name of isotope

        Returns
        -------
        isotope_dictionary : dict[ str, dict[ str, list[float] ] ]
            results for a single isotope e.g. {
            "isotope name" : {
            "activities":[], 
            "activity_uncertainties":[],
            "pathway_probabilities":[],
            "reaction_rates":[],
            "reaction_rate_uncertainty":[] }
            }
        """
        # print and save results for individual isotope activities
        # and uncerts for top 5 gamma emissions
        final_activity_list = []
        final_uncert_list = []
        intensity, energy, halflife = self._get_decay_database(isotope_name)
        for n in range(len(intensity[:5])):
            if self.json_file_data[isotope_name]['counts'][n] != 0:

                self_attenuation_factor = self._self_attenuation_correction(
                    self.json_file_data[isotope_name]['foil_material'],
                    energy[n],
                    self.json_file_data[isotope_name]['thickness_cm'],
                    self.json_file_data[isotope_name]['density_gcm3'])
                coincidence_factor = (
                    1/(self.json_file_data[isotope_name]
                       ['inv_coincidence_factor'][n]))
                uncorrected_activity = self._activity_0(
                    self.json_file_data[isotope_name]['counts'][n],
                    intensity[n],
                    energy[n],
                    measurement_distance, isotope_name, halflife)
                final_activity = (
                    coincidence_factor * self_attenuation_factor
                    * uncorrected_activity)
                counts_frac_uncert = (
                    self.json_file_data[isotope_name]['uncertainty'][n]
                    / self.json_file_data[isotope_name]['counts'][n])
                final_uncertainty = (
                    final_activity
                    * np.sqrt(counts_frac_uncert**2
                              + eff_uncert**2))
                print(f"(e={energy[n]}keV, i={intensity[n]})"
                      "activity at end of irradiation is "
                      f"{final_activity:.5e} +- {final_uncertainty:.5e} Bq")
                final_activity_list.append(final_activity)
                final_uncert_list.append(final_uncertainty)

        # print and save average reaction rates for each isotope
        # for the top gamma peak
        pathway_prob = [1]
        final_rr_list = []
        final_rr_uncert_list = []
        if "pathway_probabilities" in self.json_file_data[isotope_name]:
            pathway_prob = (self.json_file_data[isotope_name]
                            ["pathway_probabilities"])
        for p in pathway_prob:
            rr = p*self._reaction_rates(np.mean(final_activity_list),
                                        irrad_time, halflife)
            rr_uncert = p*self._reaction_rates(np.mean(final_uncert_list),
                                               irrad_time, halflife)
            final_rr_list.append(rr)
            final_rr_uncert_list.append(rr_uncert)
            print(f"Average (fraction={p}) reaction rate "
                  "over irradiation from top peak: "
                  f"{rr:.5e} +- {rr_uncert:.5e}")

        # save activites for all peaks and RRs calculated from top peak
        # as a dictionary
        isotope_dictionary = {isotope_name: {
            "activities": final_activity_list,
            "activity_uncertainties": final_uncert_list,
            "pathway_probabilities": pathway_prob,
            "reaction_rates": final_rr_list,
            "reaction_rate_uncertainty": final_rr_uncert_list
        }}
        return isotope_dictionary

    def run(self, which_isotopes, efficiency_uncert_frac, irrad_time):
        """ run analysis for all isotopes requested
        and outputs as e_results json

        Parameters
        ----------
        which_isotopes : float
            switch to control which isotopes from the data to run
        efficiency_uncert_frac : float
            Fractional uncertainty for the efficiency of the detector
        irrad_time : float
            Total irradiation time for reaction rate calculation
        """
        # set up for all isotopes requested
        open(f"{self.json_path}/e_results.json", 'w').close()
        if isinstance(which_isotopes, int):
            if which_isotopes > 0:
                isotope_run_list = list(self.json_file_data.keys()
                                        )[which_isotopes:]
            if which_isotopes < 0:
                file_isotopes = len(self.json_file_data.keys())
                isotope_run_list = list(self.json_file_data.keys()
                                        )[:file_isotopes+which_isotopes]
        else:
            if which_isotopes == 'all':
                isotope_run_list = list(self.json_file_data.keys())
            else:
                isotope_run_list = list(which_isotopes.split(" "))

        # run for all isotopes
        results_dictionary = {}
        for isotope_name in isotope_run_list:
            print(f"************ activities for {isotope_name} ************")
            measurement_distance = (self.json_file_data[isotope_name]
                                    ['detector_distance_cm'])
            eff_uncert = efficiency_uncert_frac
            results_dictionary.update(self._run_one_isotope(
                isotope_name, measurement_distance, eff_uncert, irrad_time))

        # print results as one neat json for postprocessing
        with open(f"{self.json_path}/e_results.json", 'a') as output_file:
            json.dump(results_dictionary, output_file,
                      ensure_ascii=False, indent=4)
        
