from dataclasses import dataclass
import json
from math import pi, sqrt, log, exp
import numpy as np
import actigamma as ag
from scipy.integrate import quad
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

# choose peak analysis library 
peak_library = 'root'

# folder to save the results
folder_path = f"/Users/ljb841@student.bham.ac.uk/gamma_spec/proton_hpge/experimental_activities_2025_cal/{peak_library}" 

# total irradiation time secs
irrad_time = (20+67+41)*60 + 30

# datetime timestamp for the end of your irradiation
irradiation_end = datetime(2024,3,28, 18,17,32)

# detector crystal radius in cm
detector_crystal_radius = 3.25


json_path = f"/Users/ljb841@student.bham.ac.uk/gamma_spec/proton_hpge/{peak_library}_data.json"
json_file_data = json.load(open(json_path))
specific_json_data = json_file_data['Be7']
full_isotope_list = []
for isotope in json_file_data.keys():
    full_isotope_list.append(isotope)
print(list(json_file_data.keys())[2:])

###########################################################################################
###########################################################################################
###########################################################################################
##########################################################################################

# setting library isotopes to analyse and the distance they were analysed at
if automation == 'target':
    isotope_run_list = ['Be7','Zn65']
if automation == 'foils':
    isotope_run_list = list(json_file_data.keys())[2:]
else:
    isotope_run_list = list(automation.split(" "))

#if automation in ('Be7','Zn65','target'):
#    measurement_distance = 34
#else:
#    measurement_distance = 1

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
    if isotope_name == 'Co58':
        half_life = db.gethalflife('Co58m') 
    else:
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

# solid angle approx for foils (knoll p121)
def solid_angle_disc(crystal_radius: float, distance: float, foil_radius: float) -> float:
    alpha = (foil_radius / distance)**2
    beta = (crystal_radius / distance)**2
    f_1 = ((5/16) * (beta/((1+beta)**(7/2)))) - ((35/64) * ((beta**2)/((1+beta)**(9/2))))
    f_2 = ((35/128) * (beta/((1+beta)**(9/2)))) - ((315/256) * ((beta**2)/((1+beta)**(11/2)))) + ((1155/1024) * ((beta**3)/((1+beta)**(13/2))))
    return  2 * pi * (1 - (1/((1+beta)**(1/2))) - (3/8)*((alpha*beta)/((1+beta)**(5/2))) + (alpha**2)*f_1 - (alpha**3)*f_2 )

# equation for the efficiency curves used below
def efficiency_abs(energy:float,n1:float,n2:float) -> float:
     return n1 * (energy)** (n2)

# use the measurement distance and efficiency curves to calculate activity over the live time
def activity_livetime(c,i,e) :  
    if measurement_distance == 1:
        selected_efficiency = efficiency_abs(e,8.1,-0.9209) * (solid_angle_disc(detector_crystal_radius,measurement_distance,json_file_data[isotope_name]['foil_radius_cm']) / solid_angle(detector_crystal_radius,measurement_distance))
        # louis fit with ba133,cs137,co60 (omitting ba133 81keV peak) = efficiency_abs(e,8.1,-0.9209)
        #errors=1.075,0.01871
        # fit with eu152 nov2024 = efficiency_abs(e,6.16551,-0.817)
    if measurement_distance == 6:
        selected_efficiency = efficiency_abs(e,0.6264,-0.749) # Kyle fit
    if measurement_distance == 10:
        selected_efficiency = efficiency_abs(e,0.3755,-0.765) # Kyle NEW fit
    if measurement_distance == 15:
        selected_efficiency = efficiency_abs(e,0.1966,-0.763) # Kyle fit 
    if measurement_distance == 38:
        selected_efficiency = efficiency_abs(e,0.0715,-0.8631) * (solid_angle_disc(detector_crystal_radius,measurement_distance,json_file_data[isotope_name]['foil_radius_cm']) / solid_angle(detector_crystal_radius,measurement_distance)) 
        # louis fit with ba133,cs137,co60 (omitting ba133 81keV peak) = efficiency_abs(e,0.0715,-0.8631)
        #errors=0.01523,0.03235
        # fit with eu152 nov2024 = efficiency_abs(e,0.0604,-0.778)
    activity = (c) / ((i)
        * selected_efficiency
    )
    return activity

# gamma self absorption correction factor
def self_attenuation_correction(material,E):
    xcom=np.loadtxt(f'XCOM_new/{material}.txt', skiprows=1)
    mass_coeff = np.interp(E/1000, xcom[:,0], xcom[:,1]) 
    self_att_factor = (mass_coeff * density * thickness) / (1 - exp(- mass_coeff * density * thickness))
    return self_att_factor

# integral to calculate initial activity from the measured activity over a live time
def activity_integrand(t,half_life):
    return exp(- log(2) * (t/half_life))
def activity_0(c,i,e) :
    activity = activity_livetime(c,i,e) / (quad(activity_integrand, decay_time(isotope_name), (decay_time(isotope_name)+json_file_data[isotope_name]['live_time']), args=(get_decay_database(isotope_name)[2])))
    #activity = (activity_livetime(c,i,e)/json_file_data[isotope_name]['live_time']) * exp(get_decay_database(isotope_name)[2] / (log(2) * decay_time(isotope_name)))
    return activity

# calculate reaction rate from the activity a0 at time t0 under irradiation 
def reaction_rates(a, irrad_time):
    rr_ave = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] / (1 - activity_integrand(irrad_time,get_decay_database(isotope_name)[2]))
    #rr_ave = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] / irrad_time
    rr_max = activity_0(a,get_decay_database(isotope_name)[0][0],get_decay_database(isotope_name)[1][0])[0] * (10/(8.778*128.5))
    return rr_ave,rr_max

# print and save results for individual isotope activities
for isotope_name in isotope_run_list:
    print(f"************ activities for {isotope_name} ************")
    measurement_distance = json_file_data[isotope_name]['detector_distance_cm']
    with open(f"{folder_path}/{isotope_name}_activities.txt", 'w') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:5])):
            if json_file_data[isotope_name]['counts'][n] != 0:
                print(f"(e={get_decay_database(isotope_name)[1][n]}keV, i={get_decay_database(isotope_name)[0][n]}) activity at end of irradiation is {activity_0(json_file_data[isotope_name]['counts'][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} +- {activity_0(json_file_data[isotope_name]['uncertainty'][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} Bq"
                    )
                output_file.writelines(f" e={get_decay_database(isotope_name)[1][n]} keV, i={get_decay_database(isotope_name)[0][n]}: activity at end of irradiation is {activity_0(json_file_data[isotope_name]['counts'][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} +- {activity_0(json_file_data[isotope_name]['uncertainty'][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]} Bq \n"
                                       )    

# print and save activities and uncertainties for all analysed isotopes as one nice txt
    with open(f"{folder_path}/exp_activities.txt", 'a') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:1])):
            if json_file_data[isotope_name]['counts'][n] != 0:
                output_file.writelines(f"{activity_0(json_file_data[isotope_name]['counts'][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]}\n")
    with open(f"{folder_path}/exp_uncertainties.txt", 'a') as output_file:
        for n in range(len(get_decay_database(isotope_name)[0][:1])):
            if json_file_data[isotope_name]['counts'][n] != 0:
                output_file.writelines(f"{activity_0(json_file_data[isotope_name]['uncertainty'][n],get_decay_database(isotope_name)[0][n],get_decay_database(isotope_name)[1][n])[0]}\n")

# print and save average reaction rates for each isotope in one nice txt file for unfolding
    if reaction_rate_calculator == True:
        print(f"-------- reaction rates for {isotope_name} ----------")
        pathway_prob = [1]
        if isotope_name == 'Cu64':
            pathway_prob = [0.342,0.658]
        if isotope_name == 'Cd111m':
            pathway_prob = [0.02776,0.80484,0.16737]
        if isotope_name == 'Cd115':
            pathway_prob = [0.72795,0.27205]
        if isotope_name == 'Dy157':
            pathway_prob = [0.73696,0.26304]
        for p in pathway_prob:
            print(f"Average (fraction={p}) reaction rate over irradiation:{p*reaction_rates(json_file_data[isotope_name]['counts'][0],irrad_time)[0]:.3e} +- {p*reaction_rates(json_file_data[isotope_name]['uncertainty'][0],irrad_time)[0]:.3e}")
            if automation == 'foils':
                with open(f"{folder_path}/reaction_rates.txt", 'a') as output_file:        
                    output_file.writelines(f"{p*reaction_rates(json_file_data[isotope_name]['counts'][0],irrad_time)[0]} \n")
                with open(f"{folder_path}/reaction_rate_uncertainties.txt", 'a') as output_file:    
                    output_file.writelines(f"{p*reaction_rates(json_file_data[isotope_name]['uncertainty'][0],irrad_time)[0]} \n") 
    
