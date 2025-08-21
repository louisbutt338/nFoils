"""
module for analysis of irradiation histories
look no dependencies hooray
"""

class IrradTimings:
    """ class for processing timings into fispact irradiation history
    and list of currents
    """
    def __init__(self, input_file,output_file):
        """ Initialise class
        """

        #set parameters
        self.filename = input_file
        self.output_filename = output_file

    def _parse_txt(self):
        """ Parse txt data into two arrays
        """
        input_filepath = self.filename
        with open(input_filepath,'r') as txt_data_file:
            txt_contents = txt_data_file.readlines()
            time = []
            countrate = []
            for line in txt_contents:
                time.append(float(line.split()[0]))
                countrate.append(float(line.split()[1]))
        return time,countrate

    def _proton_flux_conversion(self):
        """ return protons per second for each time interval,
        based on 12uA=190k counts
        """
        cps_to_protons_per_s = 6.24151e12 * (12/190000)
        target_protons_per_s = [
            i*cps_to_protons_per_s for i in self._parse_txt()[1] ]
        return target_protons_per_s

    def fispact_hist_writer(self):
        """ write a fispact irradiation history
        """
        with open(self.output_filename, 'w') as irrad_history_file:
            for timestep in range(len(self._parse_txt()[0])):
                irrad_history_file.writelines(
                    f"FLUX {self._proton_flux_conversion()[timestep]} \n")
                irrad_history_file.writelines(
                    "TIME 1 SECS \n")
            irrad_history_file.writelines("ATOMS")

    def run(self):
        """ run the thing how you want
        """
        #print(proton_flux_conversion())
        #fispact_hist_writer()
        #print(len(parse_txt()[0]))
        approx_current_array = [i*(12/190000) for i in self._parse_txt()[1]]
        print(approx_current_array)