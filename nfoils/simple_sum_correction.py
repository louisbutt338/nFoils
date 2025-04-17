
####################################
############ USER INPUTS ###########
####################################

# lists of energies and counts per second
#eu152_endcap = { 121.78:1541, 244.7:260, 344.28:758.2,  443.96:64.99, 778.90:172.3, 867.38:45.64, 964.08:168.6, 1085.86:121.5, 1112.07:143.3, 1408.01:178.4, 1299.14:14.29}
#eu152_far = { 121.78:16.16, 244.7:3.262, 344.28:8.735, 443.96:0.7764 , 778.90:2.189, 867.38:0.659, 964.08:2.111, 1085.86:1.425, 1112.07:1.794, 1408.01:2.331, 1299.14:0.1887}
eu152_endcap = { 121.78:2318, 244.7:257.6, 344.28:729.9,  443.96:51.36, 778.90:124.8, 867.38:27.46, 964.08:104.6, 1085.86:81.4, 1112.07:85.35, 1408.01:101.6, 1299.14:9.058}
eu152_far = { 121.78:45.63, 244.7:6.497, 344.28:13.93, 443.96:1.14 , 778.90:2.479, 867.38:0.68, 964.08:2.164, 1085.86:1.354, 1112.07:1.599, 1408.01:1.983, 1299.14:0.1928}

#am241_endcap = { 59.54:5128 }
#am241_far = { 59.54:41.24}
am241_endcap = { 59.54:17950 }
am241_far = { 59.54:350.3}

####################################
####################################

# find counts in peaks for the energy and isotope specified
def number_counts(near_dictionary,far_dictionary,energy):
    selected_near_dictionary = near_dictionary
    selected_far_dictionary = far_dictionary
    near_counts = selected_near_dictionary[energy]
    far_counts = selected_far_dictionary[energy]
    return near_counts,far_counts

# find factor for each peak at a given energy and counts in the peak. taken from simple c summing paper
def correction_factor(source_dict_near,source_dict_far,energy):
    ratio_near = number_counts(am241_endcap,am241_far,59.54)[0] / number_counts(source_dict_near,source_dict_far,energy)[0]
    ratio_far = number_counts(am241_endcap, am241_far,59.54)[1] / number_counts(source_dict_near,source_dict_far,energy)[1]
    correction_factor = ratio_near/ratio_far
    return correction_factor

correction_factors_list = []
for energy in eu152_endcap.keys():
    correction_factors_list.append(correction_factor(eu152_endcap,eu152_far,energy))
print(eu152_endcap.keys())
print(correction_factors_list)