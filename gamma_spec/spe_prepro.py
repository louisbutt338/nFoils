from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np # type: ignore
import os
import shutil

####################################
############ USER INPUTS ###########
####################################

which_foil = 'fe'

ascii_filetag = 'uBB_20s_x60_20mins'
first_file_number = 0
last_file_number = 25

input_path = '/Users/ljb841@student.bham.ac.uk/gamma_spec/deuteron_hpge/hpge_results_b03_291124/short-lived/fe_y/'
working_path = '/Users/ljb841@student.bham.ac.uk/gamma_spec/gamma_process_spectra/workshop/dli_nov24/'

# processing MAESTRO .spe file into the correct format for gamma spec code

####################################
####################################


# parses the ascii spectrum data from the selected ascii file
def parse_ascii(spectrum_number):
    filename = f"{input_path}{ascii_filetag}_{spectrum_number}.Spe"
    with open(filename,'r') as ascii_data_file:
        ascii_contents = ascii_data_file.readlines()
        ascii_header = ascii_contents[:12]
        ascii_footer = ascii_contents[8204:]
    with open(filename,'r') as ascii_data_file:
        ascii_data_strings = ascii_data_file.read().replace(" ", "").strip().split('\n')[12:8204]
        ascii_data = [int(x) for x in ascii_data_strings]
    return ascii_header,ascii_data,ascii_footer

# copies the example input file header for gamma_process_spectra and prints the desired spectrum to the new file
def spe_preprocessor(spec_numerator):
    example_path = shutil.copyfile(f"{working_path}/example_spec_format.spe", f"{working_path}/{which_foil}_data/{spec_numerator}.spe")
    with open(example_path,'r') as example_input_spectra:
        input_file = example_input_spectra.readlines()
        input_file[11-1] = f"Real Time: {int(all_ascii_headers[spec_numerator][10-1].split()[0])*(spec_numerator+1)}\n"
        input_file[12-1] = f"Live Time: {int(all_ascii_headers[spec_numerator][10-1].split()[1])*(spec_numerator+1)}\n"
        input_file[13-1] = f"Acquisition start date: {all_ascii_headers[spec_numerator][8-1].split()[0]}\n"
        input_file[14-1] = f"Acquisition start time: {all_ascii_headers[spec_numerator][8-1].split()[1]}\n"
        for j in range(8192):
            input_file[53-1+j] = f"     {j}:    {cumulative_data[j]}\n"
    with open(example_path, 'w') as example_input_spectra:
        example_input_spectra.writelines(input_file)

# automates for all the ascii files specified in the user inputs
ascii_number_array = np.arange(first_file_number,last_file_number+1)
ascii_number_array_strings = []
all_ascii_data = []
all_ascii_headers = []
for i in ascii_number_array:
    i_string = (f"{i :03d}")
    ascii_number_array_strings.append((f"{i :03d}"))
    all_ascii_headers.append(parse_ascii(f"{i :03d}")[0])
    all_ascii_data.append(parse_ascii(f"{i :03d}")[1])

    cumulative_data = [sum(x) for x in zip(*all_ascii_data)]
    print(cumulative_data[500:550])

    spe_preprocessor(i)

ascii_histogram = [sum(x) for x in zip(*all_ascii_data)]
