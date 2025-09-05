""" example for getting specific activity data from a fispact json
"""
from nfoils.fispact import JsonRetriever

# input working directory with the json file in 
upper_directory = '../analyse_fispact_json'

# input json filename (not including extension or the material bit)
json_name_format = 'example'

#input time interval in the .out file where the data u want is 
time_interval = 5

# dict of materials and the zai of the isotopes 
# include all you want to know activity for
isotope_dictionary = {
    'mn56' :{'foil':'fe','isotope_zai':250560,'isotope_state':''}
}

# loop through results - needs encapsulating
for i in isotope_dictionary.items():
    isotope = i[0]
    material = i[1]['foil']
    zai= i[1]['isotope_zai']
    isotope_state = i[1]['isotope_state']
    json_path = f'{upper_directory}/{json_name_format}.json'
    isotope_state = i[1]['isotope_state']
    retrieve_activities = JsonRetriever(
        json_path,time_interval,zai,isotope_state) 
    isotope_activity = retrieve_activities.run()
    print(f'{isotope} activity : {isotope_activity}')
