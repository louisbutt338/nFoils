from nfoils.fispact import JsonRetriever

# input working directory with the json file in 
upper_directory = '../analyse_fispact_json'

# input json filename (not including extension or the material bit)
json_name_format = 'example'

#input time interval in the .out file where the data u want is 
time_interval = 1

# dict of materials and the zai of the isotopes 
# include all you want to know activity for
isotope_dictionary = {
    'mn56' :{'foil':'fe','isotope_zai':250560,'isotope_state':''}
}

#input list of the folders with the jsons in
folders_list = ['endfb8','irdff2','tendl21']

# loop through results - needs encapsulating
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

