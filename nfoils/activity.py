"""
module for calculating activity for a whole ton of isotopes. 
thanks to Tony Turner (UKAEA) for supporting development

G Knoll, Radiation Detection and Measurement, 2010

D Arnold et al, Fundamentals of gamma spectrometry, 2018
"""
import json
from datetime import datetime
from math import pi, sqrt, log, exp
import numpy as np
import actigamma as ag
from scipy.integrate import quad
import os
import csv


class ActivityCalc:
    """ class that calculates activities for gamma spec data. 
    this gamma spec data should be written in a provided example json
    """

    def __init__(self, data_file_name, json_path, irradiation_end,
                 cal_file_name):
        """ Initialise ActivityCalc class

        Attributes
        ----------
        data_file_name : str
            name of the foil data json file
        json_path : str
            path to the directory with the json data in it
        irradiation_end : datetime object
            datetime object for the exact end of the irradiation
        cal_file_name : str
            name of the calibration data json file
        """

        # set attributes
        self.data_file_name = data_file_name
        self.json_path = json_path
        self.irradiation_end = irradiation_end
        self.cal_file_name = cal_file_name

        # load the foil data
        with open(f'{self.json_path}/{self.data_file_name}.json'
                  ) as json_datafile:
            self.json_file_data = json.load(json_datafile)

        # load the calibration data
        with open(f'{self.json_path}/{self.cal_file_name}.json'
                  ) as cal_file:
            self.cal_file_data = json.load(cal_file)

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

    def _activity_livetime(self, c, i, e, calibration_name, isotope_name):
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
        calibration_name : str
            Name of calibration in calibration_data file 
            i.e. "B03_hpge_endcap"
        isotope_name : str
            Name of isotope

        Returns
        -------
        activity : float
            the foil activity over detector livetime
        """
        # figure out the data specfic to the calibration and foil
        foil_dist = (self.cal_file_data[calibration_name]
                     ["foil_distance_cm"])
        detector_radius = (self.cal_file_data[calibration_name]
                           ["detector_radius_cm"])
        eff_values = (self.cal_file_data[calibration_name]
                      ["efficiency_equation_parameters"])
        foil_rad = self.json_file_data[isotope_name]['foil_radius_cm']

        # do the activity calculation
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

    def _activity_0(self, c, i, e, calibration_name,
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
        calibration_name : str
            Name of calibration in calibration_data file 
            i.e. "B03_hpge_endcap"
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
        activity = (self._activity_livetime(c, i, e, calibration_name,
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

    def _run_one_isotope(self, isotope_name, calibration_name, irrad_time):
        """ run analysis for one isotope

        Parameters
        ----------
        isotope_name : str
            Name of isotope
        calibration_name : str
            Name of calibration in calibration_data file 
            i.e. "B03_hpge_endcap"
        irrad_time : int
            Irradiation time in seconds

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
        # get the efficiency uncertainty for specified calibration
        eff_uncert = (self.cal_file_data[calibration_name]
                      ["efficiency_curve_fractional_uncertainty"])

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
                inv_coincidence_factor = (
                    1/(self.json_file_data[isotope_name]
                       ['coincidence_factor'][n]))
                uncorrected_activity = self._activity_0(
                    self.json_file_data[isotope_name]['counts'][n],
                    intensity[n],
                    energy[n],
                    calibration_name, isotope_name, halflife)
                final_activity = (
                    inv_coincidence_factor * self_attenuation_factor
                    * uncorrected_activity)
                counts_frac_uncert = (
                    self.json_file_data[isotope_name]['counts_uncertainty'][n]
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

    def run(self, which_isotopes, irrad_time,results_name):
        """ run analysis for all isotopes requested
        and outputs as a nice json for C/E plotting

        Parameters
        ----------
        which_isotopes : float
            switch to control which isotopes from the data to run
        irrad_time : float
            Total irradiation time for reaction rate calculation
        results_name : str
            Name of results file
        """
        # set up for all isotopes requested
        open(f"{self.json_path}/{results_name}.json", 'w').close()
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

        # get specific calibration dataset for each isotope and run
        results_dictionary = {}
        for isotope_name in isotope_run_list:
            print(f"************ activities for {isotope_name} ************")
            calibration_name = (self.json_file_data[isotope_name]
                                ["calibration"])
            results_dictionary.update(self._run_one_isotope(
                isotope_name, calibration_name, irrad_time))

        # print results as one neat json for postprocessing
        with open(f"{self.json_path}/{results_name}.json", 'a') as output_file:
            json.dump(results_dictionary, output_file,
                      ensure_ascii=False, indent=4)


class ReactionRateRetrieval:
    """ class to retrieve reaction rates for spectra-uf unfolding 
    """

    def __init__(self):
        """ Initialise RRR class

        Attributes
        ----------
        results_file_name : str
            name of the experimental results file output by ActivityCalc
        exp_folder : str
            path to the experiment directory with the results file in it
        """

    def _retrieve_rr_data(self,results_file_data):
        """ get istopes, reaction rates and uncerts for non-interfering reactions
        
        Parameters
        ----------
        results_file_data : dict
            Json dictionary object of all data in the results file 
            generated by ActivityCalc

        Returns
        -------
        isotope_list : list[str]
            list of isotopes with non-interfering reaction rates
        rr_list : list[flt]
            list of the non-interfering experimental reaction rates
        rr_u_list : list[flt]
            list of the uncertainties on the non-interfering reaction rates
        """
        # set up lists
        isotope_list = []
        rr_list = []
        rr_u_list = []

        # loop through the datasets and grab the data
        for i in results_file_data.keys():
            if len(results_file_data[i]["pathway_probabilities"]) == 1:
                isotope_list.append(
                    [i])
                rr_list.append(
                    results_file_data[i]["reaction_rates"])
                rr_u_list.append(
                    results_file_data[i]["reaction_rate_uncertainty"])
        return isotope_list, rr_list, rr_u_list
    
    def _make_new_folder(self,exp_folder):
        """ make new directory inside experimental directory to put reaction rate
        files inside

        Parameters
        ----------
        exp_folder : str
            Path to the experimental folder to work in 

        Returns
        -------
        rr_folder : str
            Path to the new folder to dump the reaction rate data into
        """
        rr_folder = os.path.join(exp_folder, "reaction_rates")
        os.makedirs(rr_folder,exist_ok=True)
        return rr_folder

    def _dump_results(self,rr_folder,isotope_list,rr_list,rr_u_list):
        """ dump results in three files for spectra-uf to read

        Parameters
        ----------
        rr_folder : str
            Path to the new folder to dump the reaction rate data into
        isotope_list : list[str]
            list of isotopes with non-interfering reaction rates
        rr_list : list[flt]
            list of the non-interfering experimental reaction rates
        rr_u_list : list[flt]
            list of the uncertainties on the non-interfering reaction rates
        """
        # list the three filenames and three datasets
        filename_list = ['isotopes','reaction_rates','reaction_rate_uncerts']
        results_list_of_lists = [isotope_list,rr_list,rr_u_list]

        # loop through the datasets and dump the data in a file for each
        for i,j in zip(filename_list,results_list_of_lists):
            open(f'{rr_folder}/{i}.csv', 'w').close()
            with open(f'{rr_folder}/{i}.csv', 'a', newline='') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerows(j)
    
    def run(self,exp_folder,results_file_name):
        """ get the reaction rate data out and dump it

        Parameters
        ----------
        exp_folder : str
            Path to the experimental folder to work in 
        results_file_name : str
            Name of the results file generated by ActivityCalc 
        """

        # load the results file 
        with open(f'{exp_folder}/{results_file_name}.json' 
                  ) as results_file:
            results_file_data = json.load(results_file)

        # retrieve data, make the folder, and dump the results inside
        isotope_list,rr_list,rr_u_list = self._retrieve_rr_data(
            results_file_data)
        rr_folder = self._make_new_folder(exp_folder)
        self._dump_results(rr_folder,isotope_list,rr_list,rr_u_list)