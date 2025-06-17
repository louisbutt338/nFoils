import numpy as np # type: ignore
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},weight='normal',size=20)

class AsciiSummer:
    def __init__(self, data_folder_path, filetag, ff_number, lf_number):
        self.folder_path = data_folder_path
        self.ascii_filetag = filetag
        self.first_file_number =  ff_number
        self.last_file_number = lf_number
    
    # parses the ascii spectrum data from the selected ascii file
    def _parse_ascii(self, spectrum_number):
        filename = f"{self.folder_path}/{self.ascii_filetag}_{spectrum_number}.Spe"
        with open(filename,'r') as ascii_data_file:
            ascii_contents = ascii_data_file.readlines()
            ascii_header = ascii_contents[:12]
            ascii_footer = ascii_contents[8204:]
        with open(filename,'r') as ascii_data_file:
            ascii_data_strings = ascii_data_file.read().replace(" ", "").strip().split('\n')[12:8204]
            ascii_data = [int(x) for x in ascii_data_strings]
        return ascii_header,ascii_data,ascii_footer

    # automates for all the ascii files specified in the user inputs
    def _loop_parser(self):
        ascii_number_array = np.arange(self.first_file_number,self.last_file_number+1)
        ascii_number_array_strings = []
        all_ascii_data = []
        for n in ascii_number_array:
            n_string = (f"{n :03d}")
            ascii_number_array_strings.append(n_string)
            all_ascii_data.append(self._parse_ascii(n_string)[1])
        ascii_histogram = [sum(x) for x in zip(*all_ascii_data)]
        return ascii_histogram,ascii_number_array_strings
    
    # plot the summed data
    def _plot_ascii(self,ascii_data):
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
        fig.savefig( 'gamma_spectrum.png', transparent=False, bbox_inches='tight')

    # writes the summed output file using the header and footer data from the FIRST ascii analysed
    def run(self):
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