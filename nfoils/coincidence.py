""" 
module for performing coincidence summing calculations and corrections
wahey no dependencies
"""

class SimpleCorrection():
    def __init__(self,source_dict_near,source_dict_far,
                 am241_dict_near,am241_dict_far):
        """ class for simple sum correction 
        with an am241 source and a multi-peak source
        """

        # set parameters
        self.source_dict_near = source_dict_near
        self.source_dict_far = source_dict_far
        self.am241_dict_near = am241_dict_near
        self.am241_dict_far = am241_dict_far

    def _number_counts(self,near_dictionary,far_dictionary,energy):
        """ find counts in peaks for the energy and isotope specified

        Parameters
        ----------
        near_dictionary : dict
            Data for the measurements close to detector 
        far_dictionary : dict
            Data for the measurements far away from detector
        """
        selected_near_dictionary = near_dictionary
        selected_far_dictionary = far_dictionary
        near_counts = selected_near_dictionary[energy]
        far_counts = selected_far_dictionary[energy]
        return near_counts,far_counts

    def _correction_factor(self,source_dict_near,source_dict_far,energy):
        """ find c.factor for a given energy/counts/peak 
        taken from simple c summing paper

        Parameters
        ----------
        near_dictionary : dict
            Data for the measurements close to detector 
        far_dictionary : dict
            Data for the measurements far away from detector
        """
        ratio_near = self._number_counts(
            self.am241_dict_near,
            self.am241_dict_far,59.54)[0] / self._number_counts(
                source_dict_near,
                source_dict_far,energy)[0]
        ratio_far = self._number_counts(
            self.am241_dict_near,
            self.am241_dict_far,59.54)[1] / self._number_counts(
                source_dict_near,
                source_dict_far,energy)[1]
        correction_factor = ratio_near/ratio_far
        return correction_factor

    def run(self):
        """ run the damn thing
        """
        correction_factors_list = []
        for energy in self.source_dict_near.keys():
            correction_factors_list.append(
                self._correction_factor(
                    self.source_dict_near,self.source_dict_far,energy))
        print(self.source_dict_near.keys())
        print(correction_factors_list)