""" 
module for performing simple coincidence summing calculations and corrections
for a multi peak source

S Kafala, Simple method for true coincidence summing correction, Journal of
    Radioanalytical and Nuclear Chemistry Articles 191 (105-114), 1995
"""

import json


class SimpleCorrection:
    """ class for simple sum correction for a multi peak source e.g. Eu152
    using a single peak source e.g. Am241
    """
    
    def __init__(self,cps_json_path):
        """ initialise SimpleCorrection class

        Attributes
        ----------
        cps_json_path : str
            path to the json data file with gamma energy (keV) and cps data
            for the multi-peak source and single-peak source at close geometry
            and at far geometry
        """

        # set attributes
        self.cps_json_name = cps_json_path

        # load the cps and energy data
        with open(f'{cps_json_path}.json'
                  ) as cps_data_file:
            self.cps_data = json.load(cps_data_file)

        # set dicts for the rest of the code
        self.multi_dict_near = self.cps_data["multi_source_near"]
        self.multi_dict_far  =  self.cps_data["multi_source_far"]
        self.single_dict_near = self.cps_data["single_source_near"]
        self.single_dict_far =  self.cps_data["single_source_far"]

    def _number_counts(self,near_dictionary,far_dictionary,energy):
        """ find counts in peaks for the energy and isotope specified

        Parameters
        ----------
        near_dictionary : dict
            Data for the measurements close to detector 
        far_dictionary : dict
            Data for the measurements far away from detector
        energy : float
            Energy value in keV for investigation

        Returns
        ----------
        near_counts : list[float]
            countrates at the near geometry, for given energies
        far_counts : list[float]
            countrates at the far geometry, for given energies
        """
        energy_str = f'{energy}'
        selected_near_dictionary = near_dictionary
        selected_far_dictionary = far_dictionary
        near_counts = selected_near_dictionary[energy_str]
        far_counts = selected_far_dictionary[energy_str]
        return near_counts,far_counts

    def _correction_factor(self,source_dict_near,source_dict_far,energy):
        """ find c.factor for a given energy/counts/peak (kafala)

        Parameters
        ----------
        source_dict_near : dict
            Data for measurements close to detector 
        source_dict_far : dict
            Data for measurements far away from detector
        energy : float
            Energy value in keV for investigation

        Returns
        ----------
        correction_factor : float
            simple correction factor for the energy specified
        """
        single_energy = [float(i) for i in self.single_dict_near.keys()]
        ratio_near = self._number_counts(
            self.single_dict_near,
            self.single_dict_far,single_energy[0])[0] / self._number_counts(
                source_dict_near,
                source_dict_far,energy)[0]
        ratio_far = self._number_counts(
            self.single_dict_near,
            self.single_dict_far,single_energy[0])[1] / self._number_counts(
                source_dict_near,
                source_dict_far,energy)[1]
        correction_factor = ratio_near/ratio_far
        return correction_factor

    def run(self):
        """ run the damn thing
        """
        correction_factors_list = []
        energy_list = []
        for energy in self.multi_dict_near.keys():
            energy_list.append(energy)
            correction_factors_list.append(
                self._correction_factor(
                    self.multi_dict_near,self.multi_dict_far,float(energy)))
        print('energies (keV):',energy_list)
        print('factors:',correction_factors_list)