import json
from math import pi, sqrt, log, exp
import numpy as np 
import actigamma as ag 
from scipy.integrate import quad 
from datetime import datetime

class ActivityCalc:
    def __init__(self, reaction_rate_calculator,automation,
            data_file_name,efficiency_uncertainty_frac,
            json_path,irrad_time,irradiation_end):
        self.reaction_rate_calculator = reaction_rate_calculator
        self.automation =  automation
        self.data_file_name = data_file_name
        self.efficiency_uncertainty_frac = efficiency_uncertainty_frac
        self.json_path = json_path
        self.irrad_time = irrad_time
        self.irradiation_end = irradiation_end
        self.json_file_data = json.load(open(
            f'{self.json_path}/{self.data_file_name}.json'))

    # function for calculating the decay time 
    # from input irradiation end and measurement start
    def _decay_time(self,isotope_name):
        dt_elements = self.json_file_data[isotope_name]['datetime']
        ts = datetime(*dt_elements)
        return (ts-self.irradiation_end).total_seconds()

    # read the pypact actigamma decay2012 database 
    # and get halflife, intensities,peaks for specified isotope
    def _get_decay_database(self,isotope_name):
        SPECTYPE = "gamma"
        db = ag.Decay2012Database()
        half_life = db.gethalflife(isotope_name) 
        intensity = db.getintensities(isotope_name,spectype=SPECTYPE)
        peak_energy_kev = 1e-3 * db.getenergies(isotope_name, spectype=SPECTYPE)
        sorted_lists = sorted(zip(peak_energy_kev, intensity),key=lambda x: x[1])
        e, i = zip(*reversed(sorted_lists))
        return(i,e,half_life)

    # solid angle approximation (knoll, p120)
    def _solid_angle(self,crystal_radius: float, distance: float) -> float:
        return 2 * pi * (1 - distance / sqrt(distance**2 + crystal_radius**2))

    # solid angle approx for foils/discs (knoll p121)
    def _solid_angle_disc(self,crystal_radius: float, 
                          distance: float, foil_radius: float) -> float:
        alpha = (foil_radius / distance)**2
        beta = (crystal_radius / distance)**2
        f_1 = ( ((5/16) * (beta/((1+beta)**(7/2)))) 
               - ((35/64) * ((beta**2)/((1+beta)**(9/2)))) )
        f_2 = ( ((35/128) * (beta/((1+beta)**(9/2)))) 
               - ((315/256) * ((beta**2)/((1+beta)**(11/2)))) 
               + ((1155/1024) * ((beta**3)/((1+beta)**(13/2)))) )
        return ( 2 * pi * (1 - (1/((1+beta)**(1/2)))
                            - (3/8)*((alpha*beta)/((1+beta)**(5/2)))
                              + (alpha**2)*f_1 - (alpha**3)*f_2 ) )

    # equation for the log-polynomial efficiency curves
    def _efficiency_abs(self,energy:float,
                        n0:float,n1:float,n2:float,n3:float) -> float:
        polynomial = ( n0 + n1*np.log(energy)**1 
                      + n2*np.log(energy)**2 + n3*np.log(energy)**3  )
        return np.exp(polynomial)

    # use the foil measurement distance and efficiency curves 
    # to calculate activity over the live time - FIX G11/B03 DISCREPANCY
    def _activity_livetime(self,c,i,e,foil_dist,isotope_name) :  
        if foil_dist == 1: #b03 hpge
            detector_radius = 3.25
            eff_values = [-23.385, 11.3457, -1.9302, 0.10197]
        if foil_dist == 0.4: #g11 bege
            detector_radius = 3.75
            eff_values = [-25.7209, 13.8415, -2.48, 0.1352]
        if foil_dist == 38: #b03 hpge
            detector_radius = 3.25
            eff_values = [-27.773, 11.13964, -1.8975, 0.10025] 
        foil_rad = self.json_file_data[isotope_name]['foil_radius_cm']
        solid_angle_ratio = ( self._solid_angle_disc(detector_radius,foil_dist,foil_rad) 
                            / self._solid_angle(detector_radius,foil_dist) )
        selected_efficiency = self._efficiency_abs(e,*eff_values) * solid_angle_ratio
        activity = c / ((i)* selected_efficiency)
        return activity

    # gamma self absorption correction factor taken from XCOM mu data
    def _self_attenuation_correction(self,material,E,thickness,density):
        xcom=np.loadtxt(f'../../data/XCOM_new/{material}.txt', skiprows=1)
        mass_coeff = np.interp(E/1000, xcom[:,0], xcom[:,1]) 
        self_att_factor = ( (mass_coeff * density * thickness) 
                           / (1 - exp(- mass_coeff * density * thickness)) )
        return self_att_factor

    # integrals for activity interpolation
    def _activity_integrand(self,t,half_life):
        return exp(- log(2) * (t/half_life))
    
    # calculate initial activity (w/o corrections) 
    # from the measured activity over a live time
    def _activity_0(self,c,i,e,measurement_distance,isotope_name):
        top_time_band = (self._decay_time(isotope_name)
                         +self.json_file_data[isotope_name]['live_time'])
        quad_integral = quad(self._activity_integrand, 
                             self._decay_time(isotope_name),top_time_band,
                             args=(self._get_decay_database(isotope_name)[2]))
        activity = (self._activity_livetime(c,i,e,measurement_distance,isotope_name)
                     / (quad_integral))[0]
        #activity = ( (activity_livetime(c,i,e)[0]
        #            /json_file_data[isotope_name]['live_time']) 
        #            * exp(get_decay_database(isotope_name)[2] 
        #                / (log(2) * decay_time(isotope_name))) )
        return activity

    # calculate reaction rate from the activity a0 at time t0 under irradiation 
    def _reaction_rates(self,a, irradiation_time,halflife):
        rr_ave = a / (1 - self._activity_integrand(irradiation_time,halflife))
        return rr_ave

    def _run_one_isotope(self,isotope_name,measurement_distance):
            # print and save results for individual isotope activities 
            # and uncerts for top 5 gamma emissions
            final_activity_list = []
            final_uncert_list = []
            intensity,energy,halflife = self._get_decay_database(isotope_name) 
            for n in range(len(self._get_decay_database(isotope_name)[0][:5])):
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
                        measurement_distance,isotope_name)
                    final_activity = (
                        coincidence_factor * self_attenuation_factor 
                        * uncorrected_activity)
                    counts_frac_uncert = (
                        self.json_file_data[isotope_name]['uncertainty'][n]
                        /self.json_file_data[isotope_name]['counts'][n])
                    final_uncertainty = (
                        final_activity 
                        * np.sqrt( counts_frac_uncert**2 
                                  + self.efficiency_uncertainty_frac**2 ))

                    print(f"(e={energy[n]}keV, i={intensity[n]}) " 
                          "activity at end of irradiation is " 
                           f"{final_activity:.5e} +- {final_uncertainty:.5e} Bq")
                    final_activity_list.append(final_activity)
                    final_uncert_list.append(final_uncertainty)

            # print and save average reaction rates for each isotope
            # for the top gamma peak
            pathway_prob = [1]
            final_rr_list = []
            final_rr_uncert_list = []
            if self.reaction_rate_calculator == True:
                if "pathway_probabilities" in self.json_file_data[isotope_name]:
                    pathway_prob = (self.json_file_data[isotope_name]
                                    ["pathway_probabilities"])
                for p in pathway_prob:
                    final_rr = p*self._reaction_rates(np.mean(final_activity_list),
                                                      self.irrad_time,halflife)
                    final_rr_uncert = p*self._reaction_rates(np.mean(final_uncert_list),
                                                             self.irrad_time,halflife)
                    final_rr_list.append(final_rr)
                    final_rr_uncert_list.append(final_rr_uncert)
                    print(f"Average (fraction={p}) reaction rate " 
                          "over irradiation from top peak: "
                           f"{final_rr:.5e} +- {final_rr_uncert:.5e}") 

            # save activites for all peaks and RRs calculated from top peak 
            # as a dictionary
            isotope_dictionary = {isotope_name: {
                "activities":final_activity_list,
                "activity_uncertainties":final_uncert_list,
                "pathway_probabilities":pathway_prob,
                "reaction_rates":final_rr_list,
                "reaction_rate_uncertainty":final_rr_uncert_list
            }}
            return isotope_dictionary
        
    def run(self):
        # set up for all isotopes requested
        open(f"{self.json_path}/e_results.json", 'w').close()
        isotope_run_list = list(self.automation.split(" "))
        if self.automation == 'target':
            isotope_run_list = list(self.json_file_data.keys())[:2]
        if self.automation == 'foils':
            isotope_run_list = list(self.json_file_data.keys())[2:]

        # run for all isotopes
        results_dictionary = {}
        for isotope_name in isotope_run_list:
            print(f"************ activities for {isotope_name} ************")
            measurement_distance = (self.json_file_data[isotope_name]
                                    ['detector_distance_cm'])
            results_dictionary.update(self._run_one_isotope(
                isotope_name,measurement_distance))

        # print results as one neat json for postprocessing
        with open(f"{self.json_path}/e_results.json", 'a') as output_file:
            json.dump(results_dictionary,output_file,ensure_ascii=False,indent=4)
