""" 
module for dealing with and processing ASCII spe files

M Gilbert, Neutron Irradiation Experiments: Automated Processing and Analysis
    of γ-spectra, Nuclear Data Sheets 119 (401-403) 2014
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

    def __init__(self, data_folder_path, filetag, file_numbers):
        """ initialise AsciiSummer class

        Attributes
        ----------
        data_folder_path : str
            path to the folder with the asciis to be summed inside
        filetag : str
            name of the ascii files WITHOUT the '_XXX' number on the end
        file_numbers : list[int]
            numbers of the first and last asciis to be summed
            e.g. [0,2] will sum ASCIIS from 000 to 002
        """

        # set attributes
        self.folder_path = data_folder_path
        self.ascii_filetag = filetag
        self.first_file_number =  file_numbers[0]
        self.last_file_number = file_numbers[1]

    def _parse_ascii(self, spectrum_number):
        """ parses the ascii spectrum data from the selected ascii file

        Parameters
        ----------
        spectrum_number : str
            number at the end of the .spe filename format as string

        Returns
        ----------
        ascii_header : list[str]
            text header of the ascii file
        ascii_data : list[int]
            data in the middle of the ascii file
        ascii_footer : list[str]
            text footer at bottom of the file
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

        Returns
        ----------
        ascii_histogram : list[int]
            summed ascii data for all the files
        ascii_number_array_strings: list[str]
            list of the numbers at the end of ascii files
            formatted as strings
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
        ascii_data : list[int]
            data in the middle of the ascii file
        """
        print("plotting summed ASCII...")
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
        filename = f"{self.folder_path}/{self.ascii_filetag}_summed.Spe"
        with open(filename,'w') as ascii_histogram_file:
            for line in self._parse_ascii(self._loop_parser()[1][0])[0]:
                ascii_histogram_file.write(line)
            for line in self._loop_parser()[0]:
                ascii_histogram_file.write(f"{line}\n")
            for line in self._parse_ascii(self._loop_parser()[1][0])[2]:
                ascii_histogram_file.write(line)

        # plot the summed ascii file in this case 
        # change 'summed' to plot another one
        summed_spe_data = self._parse_ascii('summed')[1]
        self._plot_ascii(summed_spe_data)


class AsciiPreprocessing(AsciiSummer):
    """ class for processing MAESTRO .spe file into the format for 
    UKAEA's gamma_process_spectra (gilbert)
    """

    def __init__(self, data_folder_path, filetag, ff_number, lf_number):
        """ Initialise class (inherits AsciiSummer)
        """
        super().__init__(self, data_folder_path, filetag, ff_number, lf_number)

    def _spe_preprocessor(self,spec_numerator,ascii_headers,cumulative_data):
        """ copies the example input file header for gamma_process_spectra
        and print the desired spectrum to the new format

        Parameters
        ----------
        spec_numerator : int
            numerator for the spectra you are processing
        ascii_headers : list[str]
            Headers of the ascii files
        cumulative_data : list[int]
            Summed data for the asciis
        """

        example_path = os.system(
            f"cp {self.folder_path}/{self.ascii_filetag}_000.spe "
            f"new_format_spectra/{spec_numerator}.spe")
        
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