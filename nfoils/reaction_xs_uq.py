import sandy
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},weight='normal',size=20)
import numpy as np
import csv

# user inputs
ek=sandy.energy_grids.ECCO33
library = 'jeff_33' # endfb_71 endfb_80 jendl_40u jeff_33 tendl_21 irdff_2
#material_list = [491150, 
#                 661640,
#                 791970,
#                 491150,
#                 290650,
#                 260560,
#                 130270,
#                 791970,
#                 410930,
#                 280580] 
#mt_values_list = [[102],
#                  [102],
#                  [102],
#                  [4],
#                  [103],
#                  [103],
#                  [107],
#                  [16],
#                  [16],
#                  [16]]
#reaction_labels = [r'${}^{115}$In(n,$\gamma$)',
#                r'${}^{164}$Dy(n,$\gamma$)',
#                r'${}^{197}$Au(n,$\gamma$)',
#                r"${}^{115}$In(n,n')", 
#                r'${}^{65}$Cu(n,p) *',
#                r'${}^{56}$Fe(n,p)',
#                r'${}^{27}$Al(n,$\\alpha$)', 
#                r'${}^{197}$Au(n,2n)',
#                r'${}^{93}$Nb(n,2n)',
#                r'${}^{58}$Ni(n,2n) ']

material_list = [260560]
mt_values_list = [[103]]
reaction_labels = [r'${}^{56}$Fe(n,p)']

# get covariance and standard deviation data from a material endf6 file
def _get_cov_data(material,mt_values):
    mts = mt_values
    try:
        endf_file = sandy.get_endf6_file(library, "xs", material)
        #print(endf_file)
        ekws = dict(ek=ek)
        err = endf_file.get_errorr(temperature=0,err=1,chi=False, xs=True, nubar=False, mubar=False, errorr_kws=ekws,verbose=False)["errorr33"]
    except:
        return [],np.array([])
    covariance = err.get_cov()
    std = covariance.get_std().reset_index().query("MT in @mts")
    std["MT"] = std["MT"].astype("category")
    std["STD"] *= 100
    stdev_array = np.array(std["STD"])
    print(f'-----> found reactions for material {material}: {std["MT"].cat.categories}')
    return covariance,stdev_array

# plot covariance matrix from a material endf file - NEED to CHANGE SNS plotting to matplotlib
def _plot_cov_matrix(material,mt_values):
    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
    cov =  _get_cov_data(material,mt_values)[0]
    mask = cov.data.index.get_level_values("MT") != 1
    ng = cov.data.index.get_level_values("MT")[~mask].size
    nr = cov.data.index.get_level_values("MT").unique().size - 1
    data = cov.get_corr().data.iloc[mask, mask]
    #sns.heatmap(data=data, vmin=-1, vmax=1, cmap="bwr", ax=ax)
    for i in range(1, nr):
        ax.axvline(ng * i, color="k", ls="--", lw=.5)
        ax.axhline(ng * i, color="k", ls="--", lw=.5)
    fig.tight_layout()
    fig.savefig('cov_matrix.png')

# extract stdev data from covariances and return the stdev array split by MT reactions
def _extract_stdev_data(material,mt_values):
    stdev_array = _get_cov_data(material,mt_values)[1]
    if stdev_array.size > 0:
        number_of_arrays = len(stdev_array)/(len(ek)-1)
        if number_of_arrays > 1:
            #stdev_array_split = np.array_split(stdev_array, len(stdev_array)/(len(ek)-1))
            stdev_array_split = np.split(stdev_array, number_of_arrays)
        if number_of_arrays == 1:
            stdev_array_split = [stdev_array]
        #stdev_array_transposed = stdev_array.ravel()[None]
        return stdev_array_split
    else:
        print(f'-----> reactions not found for material {material}')
        return None

# export uncert data to one csv and plot uncertainty percentages
def _export_and_plot_stdev(material_list,mt_values_list,reaction_labels):
    open('uncertainty.csv','w').close()
    np.savetxt("uncertainty_group_structure.csv", ek.ravel()[None],delimiter=',')
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(17,7),gridspec_kw={'width_ratios': [2, 3.5]},tight_layout=True)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(material_list))))

    # loop through specified materials and MT values
    for material,mts,reactions in zip(material_list,mt_values_list,reaction_labels):
        array_of_arrays = _extract_stdev_data(material,mts)
        if array_of_arrays != None:

            # plot uncertainty data
            ek_mev = [(i/1e6) for i in ek]
            for mt_iterator in range(len(array_of_arrays)):
                c=next(color)
                #print(array_of_arrays[mt_iterator])
                #print(ek_mev)
                ax1.stairs(array_of_arrays[mt_iterator], ek_mev,label=f'{reactions}',color=c)
                ax2.stairs(array_of_arrays[mt_iterator], ek_mev,label=f'{reactions}',color=c)

            # export data to one csv file
            for xs_stdev in array_of_arrays:
                #print(ek[171],xs_stdev[171])
                with open('uncertainty.csv','a',newline='') as f:
                    writer=csv.writer(f,delimiter=',' )
                    writer.writerow(xs_stdev*(1/100))
        else:
            continue
            #print(f'-----> reactions not found for material {material}')

    ax1.set_xlim(1e-8, 1e0)
    ax1.set_ylim(1e0,2e2)
    ax1.set_xscale('log')
    ax1.grid()
    ax2.set_xlim(1e0,18)
    ax2.set_ylim(1e0,2e2)
    ax2.tick_params(axis='y',left=False,labelleft=False)
    ax2.grid()  
    ax2.legend( loc="upper right", frameon=True,fontsize=18,fancybox=False,facecolor='white',framealpha=1,ncol=3)
    fig.supylabel(r"Standard deviation ($\%$)",y=0.55)
    fig.supxlabel("Neutron energy (MeV)",y=0.03)
    fig.savefig('percentage_uncert.png')

_export_and_plot_stdev(material_list,mt_values_list,reaction_labels)
