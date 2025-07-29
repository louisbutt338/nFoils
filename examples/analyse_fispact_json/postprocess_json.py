from nfoils.fispact import JsonRetriever

# input working directory with the json file in 
#upper_directory = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/model_results/approach_1_p/neutrons/calibration_model2/020725_foils_sep_ssf'
upper_directory = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/uBB/model_results/approach_1_d/neutrons/foils_seperated_010425/60mb_be7_xs'

# input json filename (not including extension or the material bit)
json_name_format = 'dli_experiment'

#input time interval in the .out file where the data u want is 
time_interval = 5531

# dict of materials and the zai of the isotopes 
# include all you want to know activity for
isotope_dictionary = {
    'mn56' :{'foil':'fe','isotope_zai':250560,'isotope_state':''},
    'au196': {'foil':'au','isotope_zai':791960,'isotope_state':''},
    'au198': {'foil':'au','isotope_zai':791980,'isotope_state':''},
    'na24': {'foil':'al','isotope_zai':110240,'isotope_state':''},
    'mg27': {'foil':'al','isotope_zai':120270,'isotope_state':''},
    'ni65': {'foil':'cu','isotope_zai':280650,'isotope_state':''},
    'cu64': {'foil':'cu','isotope_zai':290640,'isotope_state':''},
    'cd111m': {'foil':'cd','isotope_zai':481110,'isotope_state':'m'},
    'in117': {'foil':'cd','isotope_zai':491170,'isotope_state':''},
    'in115m': {'foil':'in','isotope_zai':491150,'isotope_state':'m'},
    'in116m': {'foil':'in','isotope_zai':491160,'isotope_state':'m'},
    'ni57': {'foil':'ni','isotope_zai':280570,'isotope_state':''},
    'co58': {'foil':'ni','isotope_zai':270580,'isotope_state':''},
    'co61': {'foil':'ni','isotope_zai':270610,'isotope_state':''},
    'co57': {'foil':'ni','isotope_zai':270570,'isotope_state':''},
    'dy165': {'foil':'dy','isotope_zai':661650,'isotope_state':''},
    'dy157': {'foil':'dy','isotope_zai':661570,'isotope_state':''},
    'sc44m': {'foil':'sc','isotope_zai':210440,'isotope_state':'m'},
    'y90m': {'foil':'y' ,'isotope_zai':390900,'isotope_state':'m'},
    'mo99': {'foil':'mo','isotope_zai':420990,'isotope_state':''},
    'mo99m': {'foil':'mo','isotope_zai':430990,'isotope_state':'m'},
    'mo101': {'foil':'mo','isotope_zai':431010,'isotope_state':''},
}

#input list of the folders with the jsons in
folders_list = ['endfb8','irdff2','tendl21']

for i in isotope_dictionary.items():
    for j in folders_list:
        isotope = i[0]
        material = i[1]['foil']
        lower_directory = f'{upper_directory}/{j}'
        json_path = f'{lower_directory}/{material}_{json_name_format}.json'

        zai= i[1]['isotope_zai']
        isotope_state = i[1]['isotope_state']
        retrieve_activities = JsonRetriever(
            json_path,time_interval,zai,isotope_state) 
        isotope_activity = retrieve_activities.run()
        print(f'{isotope} {j} activity : {isotope_activity}')
    print('')

