""" 
module for dealing with ASCII spe files
"""

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import rc
import os
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=20)


class AsciiSummer:
    """ class for summing time-binned ascii files together and plotting
    """
    def __init__(self, data_folder_path, filetag, ff_number, lf_number):
        """ initialise class
        """
        # set the parameters
        self.folder_path = data_folder_path
        self.ascii_filetag = filetag
        self.first_file_number =  ff_number
        self.last_file_number = lf_number

    def _parse_ascii(self, spectrum_number):
        """ parses the ascii spectrum data from the selected ascii file

        Parameters
        ----------
        spectrum_number : str
            Number at the end of the .spe filename
        """
        filename = (f"{self.folder_path}/{self.ascii_filetag}"
                    f"_{spectrum_number}.Spe")
        with open(filename,'r') as ascii_data_file:
            ascii_contents = ascii_data_file.readlines()
            ascii_header = ascii_contents[:12]
            ascii_footer = ascii_contents[8204:]
        with open(filename,'r') as ascii_data_file:
            asc_strings = ascii_data_file.read().replace(" ", "")
            asc_strings_mod = asc_strings.strip().split('\n')[12:8204]
            ascii_data = [int(x) for x in asc_strings_mod]
        return ascii_header,ascii_data,ascii_footer

    def _loop_parser(self):
        """ automate the data parsing for all the ascii files specified

        Parameters
        ----------
        spectrum_number : str
            Number at the end of the .spe filename
        """
        ascii_number_array = np.arange(self.first_file_number,
                                       self.last_file_number+1)
        ascii_number_array_strings = []
        all_ascii_data = []
        for n in ascii_number_array:
            n_string = (f"{n :03d}")
            ascii_number_array_strings.append(n_string)
            all_ascii_data.append(self._parse_ascii(n_string)[1])
        ascii_histogram = [sum(x) for x in zip(*all_ascii_data)]
        return ascii_histogram,ascii_number_array_strings
    
    def _plot_ascii(self,ascii_data):
        """ plot the summed data

        Parameters
        ----------
        ascii_data : list
            Ascii data in list format
        """
        kev_array = [i*0.41653 for i in range(len(ascii_data))]
        fig, ax1 = plt.subplots(tight_layout=True)
        ax1.set_xlabel('Gamma energy (keV)') 
        ax1.set_ylabel('Counts')
        ax1.tick_params(axis='y')
        ax1.set_xlim(0,2200)
        ax1.set_ylim(1e0,5e3)
        ax1.set_yscale("log")
        ax1.plot(kev_array, ascii_data , 'b-' )
        ax1.grid(which='major')
        fig.set_size_inches((12, 6))
        fig.savefig('gamma_spectrum.png', transparent=False,
                    bbox_inches='tight')

    def run(self):
        """ writes the summed output file 
        using the header and footer data from the FIRST ascii analysed
        """
        print('writing summed ASCII...')
        filename = f"{self.folder_path}/summed_{self.ascii_filetag}.Spe"
        with open(filename,'w') as ascii_histogram_file:
            for line in self._parse_ascii(self._loop_parser()[1][0])[0]:
                ascii_histogram_file.write(line)
            for line in self._loop_parser()[0]:
                ascii_histogram_file.write(f"{line}\n")
            for line in self._parse_ascii(self._loop_parser()[1][0])[2]:
                ascii_histogram_file.write(line)
        summed_spe_data = self._parse_ascii('000')[1]
        self._plot_ascii(summed_spe_data)


class AsciiPreprocessing(AsciiSummer):
    """ class for processing MAESTRO .spe file 
    into the format for Mark's fortran gamma spec program
    """

    def __init__(self, which_foil, ascii_filetag, ff_number, lf_number):
        """ initialise class
        """
        # set the parameters
        self.which_foil = 'fe'
        self.ascii_filetag = 'uBB_20s_x60_20mins'
        self.first_file_number = 0
        self.last_file_number = 25

    def _spe_preprocessor(self,spec_numerator,ascii_headers,cumulative_data):
        """ copies the example input file header for gamma_process_spectra
         and print the desired spectrum to the new format

        Parameters
        ----------
        spec_numerator : int
            numerator for the spectra you are processing
        ascii_headers : str?
            Headers of the ascii files
        cumulative_data : str?
            Summed data for the asciis
        """

        example_path = os.system(
            "cp example_asciis/uBB_20s_x60_20mins_000.spe "
            f"new_format_spectra/{self.which_foil}_data/{spec_numerator}.spe")
        with open(example_path,'r') as example_input_spectra:
            input_file = example_input_spectra.readlines()
            input_file[11-1] = "Real Time: "
            f"{int(ascii_headers[spec_numerator][10-1].split()[0])*(
                spec_numerator+1)}\n"
            input_file[12-1] = "Live Time: "
            f"{int(ascii_headers[spec_numerator][10-1].split()[1])*(
                spec_numerator+1)}\n"
            input_file[13-1] = "Acquisition start date: "
            f"{ascii_headers[spec_numerator][8-1].split()[0]}\n"
            input_file[14-1] = "Acquisition start time: "
            f"{ascii_headers[spec_numerator][8-1].split()[1]}\n"
            for j in range(8192):
                input_file[53-1+j] = f"     {j}:    {cumulative_data[j]}\n"
        with open(example_path, 'w') as example_input_spectra:
            example_input_spectra.writelines(input_file)

    def run(self):
        """ automates for all the ascii files specified in the user inputs
        """
        ascii_number_array = np.arange(
            self.first_file_number,self.last_file_number+1)
        ascii_number_array_strings = []
        all_ascii_data = []
        all_ascii_headers = []
        for i in ascii_number_array:
            ascii_number_array_strings.append((f"{i :03d}"))
            all_ascii_headers.append(AsciiSummer._parse_ascii(f"{i :03d}")[0])
            all_ascii_data.append(AsciiSummer._parse_ascii(f"{i :03d}")[1])
            cumulative_data = [sum(x) for x in zip(*all_ascii_data)]

            print(cumulative_data[500:550])
            self._spe_preprocessor(i,all_ascii_headers,cumulative_data)