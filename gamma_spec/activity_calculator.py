from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np # type: ignore
import actigamma as ag # type: ignore
from scipy.integrate import quad # type: ignore
from datetime import datetime

################################################################################
########### USER INPUTS ########################################################
################################################################################
################################################################################

# choose whether you want reaction rates along with the activities
reaction_rate_calculator = True

# choose whether to run all FOILS isotopes ('foils'), TARGET isotopes ('target') 
# or a specific isotope ('isotope'):
automation = 'foils'

# choose peak analysis library (root or interspec) 
peak_library = 'root'

# choose experiment (deuteron or proton)
experiment = 'proton'

# input fractional uncertainty on your eff curve fit
efficiency_uncert_frac = 0.057

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# folder for input data
json_path = f"/Users/ljb841@student.bham.ac.uk/gamma_spec/{experiment}_hpge/{peak_library}_data.json"
json_file_data = json.load(open(json_path))
# folder to save the results
folder_path = f"/Users/ljb841@student.bham.ac.uk/gamma_spec/{experiment}_hpge/experimental_activities/{peak_library}" 


# switches for the different experiments
if experiment == 'proton':
    irrad_time = (20+67+41)*60 + 30
    irradiation_end = datetime(2024,3,28, 18,17,32) 
if experiment == 'deuteron':
    irrad_time = 20*60
    irradiation_end = datetime(2024,11,29, 15,18,45)

# setting library isotopes to analyse and the distance they were analysed at
isotope_run_list = list(automation.split(" "))
if automation == 'target':
    isotope_run_list = list(json_file_data.keys())[:2]
if automation == 'foils':
    isotope_run_list = list(json_file_data.keys())[2:]

# function for calculating the decay time from input irradiation end and measurement start
def decay_time(isotope_name):
    ts = datetime(json_file_data[isotope_name]['datetime'][0],
                  json_file_data[isotope_name]['datetime'][1],
                  json_file_data[isotope_name]['datetime'][2],
                  json_file_data[isotope_name]['datetime'][3],
                  json_file_data[isotope_name]['datetime'][4],
                  json_file_data[isotope_name]['datetime'][5])
    return (ts-irradiation_end).total_seconds()


# read the pypact actigamma decay2012 database and get halflife, intensities,peaks for specified isotope
def get_decay_database(isotope_name):
    SPECTYPE = "gamma"
    db = ag.Decay2012Database()
    half_life = db.gethalflife(isotope_name) 
    intensity = db.getintensities(isotope_name,spectype=SPECTYPE)
    peak_energy_kev = (db.getenergies(isotope_name, spectype=SPECTYPE))*1e-3
    sorted_lists = sorted(zip(peak_energy_kev, intensity),key=lambda x: x[1])
    e, i = zip(*reversed(sorted_lists))
    return(i,e,half_life)

# for convenience
@dataclass
class FispactOutput:
    name: str
    activity: float

# solid angle approximation (knoll, p120)
def solid_angle(crystal_radius: float, distance: float) -> float:
    return 2 * pi * (1 - distance / sqrt(distance**2 + crystal_radius**2))

# solid angle approx for foils/discs (knoll p121)
def solid_angle_disc(crystal_radius: float, distance: float, foil_radius: float) -> float:
    alpha = (foil_radius / distance)**2
    beta = (crystal_radius / distance)**2
    f_1 = ((5/16) * (beta/((1+beta)**(7/2)))) - ((35/64) * ((beta**2)/((1+beta)**(9/2))))
    f_2 = ((35/128) * (beta/((1+beta)**(9/2)))) - ((315/256) * ((beta**2)/((1+beta)**(11/2)))) + ((1155/1024) * ((beta**3)/((1+beta)**(13/2))))
    return  2 * pi * (1 - (1/((1+beta)**(1/2))) - (3/8)*((alpha*beta)/((1+beta)**(5/2))) + (alpha**2)*f_1 - (alpha**3)*f_2 )

# equation for the log-polynomial efficiency curves
def efficiency_abs(energy:float,n0:float,n1:float,n2:float,n3:float) -> float:
    polynomial = n0 + n1*np.log(energy)**1 + n2*np.log(energy)**2 + n3*np.log(energy)**3  
    return np.exp(polynomial)

# use the measurement distance and efficiency curves to calculate activity over the live time
def activity_livetime(c,i,e) :  
    if measurement_distance == 1: #b03 hpge model
        detector_crystal_radius = 3.25
        eff_values = [-23.491, 11.4696, -1.952, 0.103] #exp
        #eff_values = [-15.197, 7.983, -1.456, 0.0789] #model
    if measurement_distance == 0.5: #g11 bege model SOON
        detector_crystal_radius = 3.75
        eff_values = [0.0, 4e-16, -2e-12, 4e-9 ] # CHANGE ASAP WHEN G11 MODEL FIXED
    if measurement_distance == 38: #b03 hpge model
        detector_crystal_radius = 3.25
        eff_values = [-27.257, 10.917, -1.8563, 0.09776] #exp
        #eff_values = [-23.169, 8.828, -1.491, 0.076] #model
    solid_angle_ratio = solid_angle_disc(detector_crystal_radius,measurement_distance,json_file_data[isotope_name]['foil_radius_cm']) / solid_angle(detector_crystal_radius,measurement_distance)
    selected_efficiency = efficiency_abs(e,eff_values[0],eff_values[1],eff_values[2],eff_values[3]) * solid_angle_ratio
    activity = (c) / ((i)
        * selected_efficiency
    )
    return activity

# gamma self absorption correction factor taken from XCOM mu data
def self_attenuation_correction(material,E,thickness,density):
    xcom=np.loadtxt(f'XCOM_new/{material}.txt', skiprows=1)
    mass_coeff = np.interp(E/1000, xcom[:,0], xcom[:,1]) 
    self_att_factor = (mass_coeff * density * thickness) / (1 - exp(- mass_coeff * density * thickness))
    return self_att_factor

# integrals to calculate initial activity (w/o corrections) from the measured activity over a live time
def activity_integrand(t,half_life):
    return exp(- log(2) * (t/half_life))
def activity_0(c,i,e) :
    quad_integral = quad(activity_integrand, decay_time(isotope_name), (decay_time(isotope_name)+json_file_data[isotope_name]['live_time']), args=(get_decay_database(isotope_name)[2]))
    activity = (activity_livetime(c,i,e) / (quad_integral))[0]
    #activity = (activity_livetime(c,i,e)[0]/json_file_data[isotope_name]['live_time']) * exp(get_decay_database(isotope_name)[2] / (log(2) * decay_time(isotope_name)))
    return activity

# calculate reaction rate from the activity a0 at time t0 under irradiation 
def reaction_rates(a, irrad_time):
    rr_ave = a / (1 - activity_integrand(irrad_time,get_decay_database(isotope_name)[2]))
    return rr_ave

# run for all isotopes requested
open(f"{folder_path}/exp_activities.txt", 'w').close()
open(f"{folder_path}/exp_uncertainties.txt", 'w').close()
open(f"{folder_path}/reaction_rates.txt", 'w').close()
open(f"{folder_path}/reaction_rate_uncertainties.txt", 'w').close()
for isotope_name in isotope_run_list:
    print( "*********************************************")
    print(f"************ activities for {isotope_name} ************")
    print( "*********************************************")
    measurement_distance = json_file_data[isotope_name]['detector_distance_cm']

# print and save results for individual isotope activities and uncerts for top 5 gamma emissions
    final_activity_list = []
    final_uncert_list = []
    with open(f"{folder_path}/{isotope_name}_activities.txt", 'w') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:5])):
            if json_file_data[isotope_name]['counts'][n] != 0:

                self_attenuation_factor = self_attenuation_correction(json_file_data[isotope_name]['foil_material'],get_decay_database(isotope_name)[1][n],json_file_data[isotope_name]['thickness_cm'],json_file_data[isotope_name]['density_gcm3'])
                coincidence_factor = 1/(json_file_data[isotope_name]['inv_coincidence_factor'][n])
                uncorrected_activity = activity_0(json_file_data[isotope_name]['counts'][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])
                final_activity = coincidence_factor * self_attenuation_factor * uncorrected_activity
                counts_frac_uncertainty = json_file_data[isotope_name]['uncertainty'][n]/json_file_data[isotope_name]['counts'][n]
                final_uncertainty = final_activity *  np.sqrt( counts_frac_uncertainty**2 + efficiency_uncert_frac**2 )

                print(f"(e={get_decay_database(isotope_name)[1][n]}keV, i={get_decay_database(isotope_name)[0][n]}) activity at end of irradiation is {final_activity:.5e} +- {final_uncertainty:.5e} Bq")
                output_file.writelines(f" e={get_decay_database(isotope_name)[1][n]} keV, i={get_decay_database(isotope_name)[0][n]}: activity at end of irradiation is {final_activity} +- {final_uncertainty} Bq \n")
                final_activity_list.append(final_activity)
                final_uncert_list.append(final_uncertainty)

# save activities and uncertainties FOR THE TOP GAMMA EMISSION for all analysed isotopes as one nice txt
    with open(f"{folder_path}/exp_activities.txt", 'a') as output_file:
        output_file.writelines(f"{final_activity_list[0]}\n")
    with open(f"{folder_path}/exp_uncertainties.txt", 'a') as output_file:
        output_file.writelines(f"{final_uncert_list[0]}\n")

# print and save average reaction rates for each isotope in one nice txt file for unfolding
    if reaction_rate_calculator == True:
        pathway_prob = [1]
        if isotope_name == 'Cu64':
            pathway_prob = [0.342,0.658]
        if isotope_name == 'Cd111m':
            pathway_prob = [0.02776,0.80484,0.16737]
        if isotope_name == 'Cd115':
            pathway_prob = [0.72795,0.27205]
        if isotope_name == 'Dy157':
            pathway_prob = [0.73696,0.26304]
        if isotope_name == 'Co61':
            pathway_prob == [0.78093,0.20468,0.0138]
        for p in pathway_prob:
            final_rr = p*reaction_rates(final_activity_list[0],irrad_time)
            final_rr_uncert = p*reaction_rates(final_uncert_list[0],irrad_time)
            print(f"Average (fraction={p}) reaction rate over irradiation from top peak: {final_rr:.5e} +- {final_rr_uncert:.5e}")
            if automation == 'foils':
                with open(f"{folder_path}/reaction_rates.txt", 'a') as output_file:        
                    output_file.writelines(f"{final_rr} \n")
                with open(f"{folder_path}/reaction_rate_uncertainties.txt", 'a') as output_file:    
                    output_file.writelines(f"{final_rr_uncert} \n") 
    
