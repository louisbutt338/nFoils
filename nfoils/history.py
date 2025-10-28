"""
module for analysis/postpro of irradiation histories
look no dependencies hooray
"""


class RateConversion:
    """ class for processing rate data timings txt into converted units
    """

    def __init__(self, input_file):
        """ Initialise IrradTimings class

        Attributes
        ----------
        input_file : str
            path to the txt data file with time and countrate as two columns
        """

        #set attributes
        self.filename = input_file

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

    def get_new_data(self,rate_to_new):
        """ get new data and timings from the rate data 

        Parameters
        ----------
        rate_to_new : float
            float value for convert from txt rate data 
            directly to new units

        Returns
        -------
        timing_list : list[float]
            List of time widths in s
        current_list : list[float]
            List of converted data from the rate data
        """
        print('calculating converted data...')
        timing_list = self._parse_txt()[0]
        new_data_list = [i*rate_to_new for i in self._parse_txt()[1]]
        return timing_list,new_data_list


class IrradTimeline:
    """ class for generating irradiation timelines for charged particles
    on targets. takes in list of currents and timings
    """

    def __init__(self, current_list,time_list):
        """ Initialise IrradTimeline class

        Attributes
        ----------
        current_list : list[float]
            list of currents in uA
        time_list : list[float]
            list of timings in s
        """

        #set attributes
        self.current_list = current_list
        self.time_list = time_list

    def _flux_conversion(self,pps_per_ua):
        """ return charged particles per second for each time interval,
        using an input conversion

        Parameters
        ----------
        pps_per_ua : float
            float value for converting from currents in uA to
            charged particles per s

        Returns
        -------
        target_particles_per_s : list[float]
            List of particles per s for each time interval
        """
        target_particles_per_s = [
            i*pps_per_ua for i in self.current_list]
        return target_particles_per_s

    def fispact_hist_writer(self,pps_per_ua,output_file):
        """ write a fispact irradiation history to file, in units of 
        charged particles per s

        Parameters
        ----------
        pps_per_ua : float
            float value for converting from currents in uA to
            charged particles per s
        output_file : str
            desired filepath to the output fispact irradiation history
        """
        print('writing fispact irradiation history...')
        output_txt_file = f"{output_file}.txt"
        with open(output_txt_file, 'w') as irrad_history_file:
            for i in range(len(self.current_list)):
                irrad_history_file.writelines(
                    f"FLUX {self._flux_conversion(pps_per_ua)[i]} \n")
                irrad_history_file.writelines(
                    f"TIME {self.time_list[i]} SECS \n")
            irrad_history_file.writelines("ATOMS")
