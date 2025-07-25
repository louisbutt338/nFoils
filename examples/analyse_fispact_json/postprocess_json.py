from nfoils.fispact_tools import JsonRetriever

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
    'fe':{'isotope_zai':[250560]},
    'au':{'isotope_zai':[790196]},
    'au':{'isotope_zai':[790198]},
    'al':{'isotope_zai':[110240]},
    'al':{'isotope_zai':[al270]},
    'cu':{'isotope_zai':[ni650]},
    'cu':{'isotope_zai':[cu640]},
    'cd':{'isotope_zai':[cd111],'isotope_state':'m'},
    'cd':{'isotope_zai':[in117]},
    'in':{'isotope_zai':[in115],'isotope_state':'m'},
    'in':{'isotope_zai':[in116],'isotope_state':'m'},
    'ni':{'isotope_zai':[ni570]},
    'ni':{'isotope_zai':[co580]},
    'ni':{'isotope_zai':[co610]},
    'ni':{'isotope_zai':[co570]},
    'dy':{'isotope_zai':[dy165]},
    'dy':{'isotope_zai':[dy157]},
    'sc':{'isotope_zai':[sc440],'isotope_state':'m'},
    'y':{'isotope_zai':[y90],'isotope_state':'m'},
    'mo':{'isotope_zai':[mo99]},
    'mo':{'isotope_zai':[tc99],'isotope_state':'m'},
    'mo':{'isotope_zai':[tc101]},

}

#input list of the folders with the jsons in
folders_list = ['endfb8','irdff2','tendl21']

for i in isotope_dictionary.keys():
    for j in folders_list:
        lower_directory = f'{upper_directory}/{j}'
        json_path = f'{lower_directory}/{i}_{json_name_format}.json'

        for k in isotope_dictionary[i]['isotope_zai']:
            if not isotope_dictionary[i]['isotope_state']:
                isotope_state = ''
            else:
                isotope_state = isotope_dictionary[i]['state']
            retrieve_activities = JsonRetriever(
                json_path,time_interval,k,isotope_state) 
            isotope_activity = retrieve_activities.run()
            print(f'({i} zai {k} activity) "{j}_values" : [{isotope_activity},0],')
    print('')

