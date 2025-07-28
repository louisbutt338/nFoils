import sandy
import os
import json
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},weight='normal',size=20)
import numpy as np
import csv

os.environ['NJOY'] = '/Users/ljb841@student.bham.ac.uk/NJOY2016/bin/njoy'

# user inputs
ek=sandy.energy_grids.VITAMINJ175
library = 'tendl_21' # endfb_71 endfb_80 jendl_40u jeff_33  tendl_21
data_file_name = 'foil_data_wo_in'
reaction_labels = [#r'${}^{115}$In(n,$\gamma$)',
                r'${}^{164}$Dy(n,$\gamma$)',
                r'${}^{197}$Au(n,$\gamma$)',
                r"${}^{115}$In(n,n')", 
                r'${}^{65}$Cu(n,p) *',
                r'${}^{56}$Fe(n,p)',
                r'${}^{27}$Al(n,$\alpha$)', 
                r'${}^{197}$Au(n,2n)',
                r'${}^{93}$Nb(n,2n)',
                r'${}^{58}$Ni(n,2n) ']

#reaction_labels = [r'${}^{56}$Fe(n,p)']
#reaction_labels = [r'${}^{115}$In(n,$\gamma$)']

# run the sandy get_endf routine
def _get_endf_file(material):
    endf_file = sandy.get_endf6_file(library, "xs", material)
    return endf_file

# get covariance, standard deviation, from a material endf6 file
def _get_errorr_data(material,mt_value):
    mt = [mt_value]
    try:
        # use the get_errorr function to grab cov,std data
        endf_file = _get_endf_file(material)
        ekws = dict(ek=ek)
        err = endf_file.get_errorr(temperature=0,err=1,chi=False, nubar=False, prod=False,mubar=False, errorr_kws=ekws,verbose=False)["errorr33"]
    except:
        print(f'-----> reactions not found for MAT {material}: MT {mt_value}')
        return [],np.array([])
    covariance = err.get_cov()
    std = covariance.get_std().reset_index().query("MT in @mt")
    xs = err.get_xs(mt=mt).data.to_numpy()
    std["MT"] = std["MT"].astype("category")
    std["STD"] *= 100
    stdev_array = np.array(std["STD"])
    xs_array = np.array([j for i in xs for j in i])
    print(f'-----> found reactions for MAT {material}: MT {mt_value}')
    return covariance,stdev_array,xs_array

def _get_gendf_data(material,mt_value):
    mt = [mt_value]
    try:
        # use the get_gendf function to grab xs data
        endf_file = _get_endf_file(material)
        ekws = dict(ek=ek,nuclide_production=True,iwt=3)
        gendf = endf_file.get_gendf(minimal_processing=True, err=0.005, temperature=0,groupr_kws=ekws)
    except:
        print(f'-----> reactions not found for MAT {material}: MT {mt_value}')
        return [],np.array([])
    xs = gendf.get_xs(mt=mt).data.to_numpy()
    xs_array = np.array([j for i in xs for j in i])
    print(f'-----> found reactions for MAT {material}: MT {mt_value}')
    xs_final = xs_array
    return xs_final


# extract data from the top function and return the stdev or xs array split by MT reactions
def _extract_array_data(data_array):
    if data_array.size > 0:
        number_of_arrays = len(data_array)/(len(ek)-1)
        if number_of_arrays > 1:
            #data_array_split = np.array_split(data_array, len(data_array)/(len(ek)-1))
            data_array_split = np.split(data_array, number_of_arrays)
        if number_of_arrays == 1:
            data_array_split = [data_array]
        #data_array_transposed = data_array.ravel()[None]
        return data_array_split
    
def _calculate_response_function(cross_section,density, mass,abundance,atomic_mass,thickness):
    foil_volume = mass/density
    atom_density = (abundance * density * 1e-24 * 6.022e23)/atomic_mass
    ss_correction_factor =  1 #( (1-np.exp(-atom_density*cross_section*thickness))/ (atom_density*cross_section*thickness) )
    response_function = atom_density*foil_volume*cross_section *ss_correction_factor
    response_function[np.isnan(response_function)] = 0
    return response_function

# export uncert data to one csv and plot uncertainty percentages
def _export_and_plot_stdev(material_list,mt_values_list,reaction_labels):
    open('uncertainty.csv','w').close()
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(18,8),
                                  gridspec_kw={'width_ratios': [2, 3.5]},tight_layout=True)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(material_list))))
    # loop through specified materials and MT values
    for material,mt,reaction in zip(material_list,mt_values_list,reaction_labels):
        nuclear_data = _get_errorr_data(material,mt)
        array_of_arrays = _extract_array_data(nuclear_data[1])
        if array_of_arrays != None:
            # plot uncertainty data
            ek_mev = [(i/1e6) for i in ek]
            for mt_iterator in range(len(array_of_arrays)):
                c=next(color)
                ax1.stairs(array_of_arrays[mt_iterator], ek_mev,label=f'{reaction}',color=c,lw=1.5)
                ax2.stairs(array_of_arrays[mt_iterator], ek_mev,label=f'{reaction}',color=c,lw=1.5)
            # export data to one csv file
            for xs_stdev in array_of_arrays:
                with open('uncertainty.csv','a',newline='') as f:
                    writer=csv.writer(f,delimiter=',' )
                    writer.writerow(xs_stdev*(1/100))
        else:
            continue
    ax1.set_xlim(1e-8, 1e0)
    ax1.set_ylim(1e0,2e2)
    ax1.set_xscale('log')
    ax1.grid()
    ax2.set_xlim(1e0,18)
    ax2.set_ylim(1e0,2e2)
    ax2.tick_params(axis='y',left=False,labelleft=False)
    ax2.grid()  
    ax2.legend( loc="upper right", frameon=True,fontsize=18,fancybox=False,
               facecolor='white',framealpha=1,ncol=3)
    fig.supylabel(r"Standard deviation ($\%$)",y=0.55)
    fig.supxlabel("Neutron energy (MeV)",y=0.03)
    fig.savefig('percentage_uncert.png')

# export xs data to one csv and plot 
def _export_and_plot_xs(
        material_list,mt_list,density_list,
        mass_list,abundance_list,
        atomic_mass_list,labels_list,thickness_list):
    open('response_function.csv','w').close()
    np.savetxt("group_structure.csv", ek.ravel()[None],delimiter=',')
    fig, (ax1,ax2) = plt.subplots(
        1,2,figsize=(18,8),
        gridspec_kw={'width_ratios': [2, 3.5]},
        tight_layout=True)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(material_list))))
    # loop through specified materials and MT values
    for (material,mt,density,mass,abundance,
         atomic_mass,reaction_label,thickness) in zip(
             material_list,mt_list,density_list,
             mass_list,abundance_list,atomic_mass_list,
             labels_list,thickness_list):
        nuclear_data = _get_gendf_data(material,mt)
        array_of_arrays = _extract_array_data(nuclear_data)

        # reaction dependent corrections
        if material==491150 and mt==4:
            array_of_arrays = [(i/8) for i in array_of_arrays]
        if material ==410930 and mt==16:
            array_of_arrays = [(i/3) for i in array_of_arrays]

        if array_of_arrays != None:
            # plot data
            ek_mev = [(i/1e6) for i in ek]
            for mt_iterator in range(len(array_of_arrays)):
                c=next(color)
                cross_section = array_of_arrays[mt_iterator]
                response_function = _calculate_response_function(
                    cross_section,density, mass,abundance,
                    atomic_mass,thickness)
                ax1.stairs(response_function, ek_mev,label=f'{reaction_label}',color=c,lw=1.5)
                ax2.stairs(response_function, ek_mev,label=f'{reaction_label}',color=c,lw=1.5)
            # export data to one csv file
            #for xs in response_function:
            with open('response_function.csv','a',newline='') as f:
                writer=csv.writer(f,delimiter=',' )
                writer.writerow(response_function)
        else:
            continue
    ax1.set_xlim(1e-8, 1e0)
    ax1.set_ylim(1e-11,1e3)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.grid()
    ax2.set_xlim(1e0,18)
    ax2.set_ylim(1e-11,1e3)
    ax2.set_yscale('log')
    ax2.tick_params(axis='y',left=False,labelleft=False)
    ax2.grid()  
    ax2.legend( loc="upper right", frameon=True,fontsize=18,fancybox=False,
               facecolor='white',framealpha=1,ncol=3)
    fig.supylabel("Response function Rn(E) (cm$^2$)",y=0.55)
    fig.supxlabel("Neutron energy (MeV)",y=0.03)
    fig.savefig('response_function.png')

def run():
    json_file_data = json.load(open(f'{data_file_name}.json'))
    material_list =   [x['mat_number'] for x in json_file_data.values()]
    mt_list =         [x['mt_value'] for x in json_file_data.values() ]
    density_list =    [x['density_gcm3'] for x in json_file_data.values()]
    mass_list =       [x['mass_g'] for x in json_file_data.values()]
    abundance_list =  [x['isotope_abundance'] for x in json_file_data.values()]
    atomic_mass_list= [x['foil_atomic_mass'] for x in json_file_data.values()]
    thickness_list =  [x['thickness_cm'] for x in json_file_data.values()]
    labels_list = reaction_labels
    _export_and_plot_stdev(material_list,mt_list,labels_list)
    #_export_and_plot_xs(
    #    material_list,mt_list,density_list,
    #    mass_list,abundance_list,atomic_mass_list,
    #    labels_list,thickness_list)
run()
