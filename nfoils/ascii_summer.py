import numpy as np # type: ignore

class AsciiSummer:
    def __init__(self, data_folder_path, filetag, ff_number, lf_number):
        self.folder_path = data_folder_path
        self.ascii_filetag = filetag
        self.first_file_number =  ff_number
        self.last_file_number = lf_number
    
    # parses the ascii spectrum data from the selected ascii file
    def parse_ascii(self, spectrum_number):
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
    def loop_parser(self):
        ascii_number_array = np.arange(self.first_file_number,self.last_file_number+1)
        ascii_number_array_strings = []
        all_ascii_data = []
        for n in ascii_number_array:
            n_string = (f"{n :03d}")
            ascii_number_array_strings.append(n_string)
            all_ascii_data.append(self.parse_ascii(n_string)[1])
        ascii_histogram = [sum(x) for x in zip(*all_ascii_data)]
        return ascii_histogram,ascii_number_array_strings

    # writes the summed output file using the header and footer data from the FIRST ascii analysed
    def write_ascii(self):
        print('writing summed ASCII...')
        filename = f"{self.folder_path}/TEST_summed_{self.ascii_filetag}.Spe"
        with open(filename,'w') as ascii_histogram_file:
            for line in self.parse_ascii(self.loop_parser()[1][0])[0]:
                ascii_histogram_file.write(line)
            for line in self.loop_parser()[0]:
                ascii_histogram_file.write(f"{line}\n")
            for line in self.parse_ascii(self.loop_parser()[1][0])[2]:
                ascii_histogram_file.write(line)