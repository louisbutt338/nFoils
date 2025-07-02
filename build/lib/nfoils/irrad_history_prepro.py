
working_dir = "/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/diamond_prepro"
filename = "Diamond_rate.txt"
output_filename = "irad_history_dli.txt"

# acquisition start time was 29/11/24 13:46:35
# txt data is time(s) and count rate (1/s)

# parse txt data into two arrays
def parse_txt():
    input_filepath = f"{working_dir}/{filename}"
    with open(input_filepath,'r') as txt_data_file:
        txt_contents = txt_data_file.readlines()
        time = []
        countrate = []
        for line in txt_contents:
            time.append(float(line.split()[0]))
            countrate.append(float(line.split()[1]))
    return time,countrate

# return protons per second for each time interval, based on 12uA=190k counts
def proton_flux_conversion():
    cps_to_protons_per_s = 6.24151e12 * (12/190000)
    target_protons_per_s = [i*cps_to_protons_per_s for i in parse_txt()[1] ]
    return target_protons_per_s

#print(proton_flux_conversion())

def fispact_hist_writer():
    output_filepath = f"{working_dir}/{output_filename}"
    with open(output_filepath, 'w') as irrad_history_file:
        for timestep in range(len(parse_txt()[0])):
            irrad_history_file.writelines(f"FLUX {proton_flux_conversion()[timestep]} \n")
            irrad_history_file.writelines("TIME 1 SECS \n")
        irrad_history_file.writelines("ATOMS")

#fispact_hist_writer()

approx_current_array = [i*(12/190000) for i in parse_txt()[1] ]
print(len(parse_txt()[0]))