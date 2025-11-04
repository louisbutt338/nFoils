""" get specific activity data for specified isotopes from a fispact json
 at a specified time interval. 
 then simulate a gamma spectrum for each isotope and activity
"""
from nfoils.fispact import JsonRetriever,GammaSpectrumModel

# path to a json generated during a fispact simulation
json_path = 'example'

# initialise the two classes
retrieve_json = JsonRetriever(json_path) 
gamma_spec = GammaSpectrumModel()

# input fispact irradiation history time interval that you want 
# to interrogate. May be useful to check .in and .out files
time_interval = 8

# input the isotope name and zai number of all the isotopes you want
# to analyse, in dictionary form
isotope_dictionary = {
    'mn56' :{'isotope_zai':250560}
}

# get activities for all the isotopes requested
for i in isotope_dictionary.items():
    zai = i[1]['isotope_zai']
    isotope_act = retrieve_json.get_acts(time_interval,zai)
    
    # print activities for all the isotopes
    isotope = i[0]
    print(f'{isotope} activity : {isotope_act} Bq')

    # plot a basic simulated gamma spectrum for each isotope
    isotope_ag_form = isotope.title()
    gamma_spec.plot_gamma(isotope_ag_form,isotope_act)
