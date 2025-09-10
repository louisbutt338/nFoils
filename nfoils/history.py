"""
module for analysis/postpro of irradiation histories
look no dependencies hooray
"""


class IrradTimings:
    """ class for processing timings into fispact irradiation history
    and list of currents
    """

    def __init__(self, input_file,output_file):
        """ Initialise IrradTimings class

        Attributes
        ----------
        input_file : str
            path to the txt data file with time and countrate as two columns
        output_file : str
            output filename to dump a fispact irradiation history into
        """

        #set attributes
        self.filename = input_file
        self.output_filename = output_file

    def _parse_txt(self):
        """ Parse txt data into time and countrate lists

        Returns
        -------
        time : list[float]
            List of times
        countrate : list[float]
            List of counts
        """
        input_filepath = f"{self.filename}.txt"
        with open(input_filepath,'r') as txt_data_file:
            txt_contents = txt_data_file.readlines()
            time = []
            countrate = []
            for line in txt_contents:
                time.append(float(line.split()[0]))
                countrate.append(float(line.split()[1]))
        return time,countrate

    def _flux_conversion(self,cps_to_particles_per_s):
        """ return charged particles per second for each time interval,
        using an input conversion

        Parameters
        ----------
        cps_to_particles_per_s : float
            float value for converting from cps to particles per s

        Returns
        -------
        target_particles_per_s : list[float]
            List of protons per s for each time interval
        """
        target_particles_per_s = [
            i*cps_to_particles_per_s for i in self._parse_txt()[1] ]
        return target_particles_per_s

    def fispact_hist_writer(self,cps_to_pps):
        """ write a fispact irradiation history to file

        Parameters
        ----------
        cps_to_pps : float
            float value for converting from rate data cps to
            charged particles per s
        """
        print('writing fispact irradiation history...')
        output_txt_file = f"{self.output_filename}.txt"
        with open(output_txt_file, 'w') as irrad_history_file:
            for timestep in range(len(self._parse_txt()[0])):
                irrad_history_file.writelines(
                    f"FLUX {self._flux_conversion(cps_to_pps)[timestep]} \n")
                irrad_history_file.writelines(
                    "TIME 1 SECS \n")
            irrad_history_file.writelines("ATOMS")

    def current_printer(self,rate_to_current):
        """ print currents from the rate data 

        Parameters
        ----------
        rate_to_current : float
            float value for convert from txt rate data 
            directly to charged particle current
        """
        print('getting list of converted currents...')
        current_array = [i*rate_to_current for i in self._parse_txt()[1]]
        print(current_array)