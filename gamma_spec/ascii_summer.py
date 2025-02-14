from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np # type: ignore


####################################
############ USER INPUTS ###########
####################################


ascii_filetag = 'ubb_400s_x60_6dot7hrs'
first_file_number = 0
last_file_number = 55

folder_path = '/Users/ljb841@student.bham.ac.uk/gamma_spec/deuteron_hpge/hpge_results_g11_291124/long-lived/y/'

# writes a summed ASCII file from the input ASCII files in the array. Takes the header and footer parameters 
# (i.e. livetime, timings, calibration) from the first file in the array, so this will need to be edited 

####################################
####################################

# parses the ascii spectrum data from the selected ascii file
def parse_ascii(spectrum_number):
    filename = f"{folder_path}{ascii_filetag}_{spectrum_number}.Spe"
    with open(filename,'r') as ascii_data_file:
        ascii_contents = ascii_data_file.readlines()
        ascii_header = ascii_contents[:12]
        ascii_footer = ascii_contents[8204:]
    with open(filename,'r') as ascii_data_file:
        ascii_data_strings = ascii_data_file.read().replace(" ", "").strip().split('\n')[12:8204]
        ascii_data = [int(x) for x in ascii_data_strings]
    return ascii_header,ascii_data,ascii_footer

# automates for all the ascii files specified in the user inputs
ascii_number_array = np.arange(first_file_number,last_file_number+1)
ascii_number_array_strings = []
all_ascii_data = []
for n in ascii_number_array:
    n_string = (f"{n :03d}")
    ascii_number_array_strings.append(n_string)
    all_ascii_data.append(parse_ascii(n_string)[1])
ascii_histogram = [sum(x) for x in zip(*all_ascii_data)]

# writes the summed output file using the header and footer data from the FIRST ascii analysed
def write_ascii():
    print('writing summed ASCII...')
    filename = f"{folder_path}summed_{ascii_filetag}.Spe"
    with open(filename,'w') as ascii_histogram_file:
        for line in parse_ascii(ascii_number_array_strings[0])[0]:
            ascii_histogram_file.write(line)
        for line in ascii_histogram:
            ascii_histogram_file.write(f"{line}\n")
        for line in parse_ascii(ascii_number_array_strings[0])[2]:
            ascii_histogram_file.write(line)
    
write_ascii()