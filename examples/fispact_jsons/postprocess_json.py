""" get specific activity data for specified isotopes from a fispact json
 at a specified time interval. simulate a gamma spectrum for each isotope
"""
from nfoils.fispact import JsonRetriever,GammaSpectrumModel

# path to the fispact generated json
json_path = 'example'

# initialise objects
retrieve_json = JsonRetriever(json_path) 
gamma_spec = GammaSpectrumModel()

#input time interval in the .out file where the data u want is 
time_interval = 8

# dict of zai of the isotopes
# include all you want to know activity for
isotope_dictionary = {
    'mn56' :{'isotope_zai':250560}
}

# get all requested activities
for i in isotope_dictionary.items():
    zai = i[1]['isotope_zai']
    isotope_act = retrieve_json.get_acts(time_interval,zai)
    
    # print activity results
    isotope = i[0]
    print(f'{isotope} activity : {isotope_act} Bq')

    # plot a basic simulated gamma spectrum
    isotope_ag_form = isotope.title()
    gamma_spec.plot_gamma(isotope_ag_form,isotope_act)
