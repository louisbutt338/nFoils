""" example for plotting spectrum after unfolding
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("font", **{"family":"sans-serif", "sans-serif":["Helvetica"]},
   weight='normal',size=20)
from bfoils.unfold import BayesianUnfolding

# get mean spectra and max/min spectra
spec1,uncert1 = np.loadtxt('4p_175_spectrum.txt',
                       unpack=True)
spec1max = [(i+j) for i,j in zip(spec1,uncert1)]
spec1min = [(i-j) for i,j in zip(spec1,uncert1)]
spec2,uncert2 = np.loadtxt('gravel_175_spectrum.txt',delimiter=',')
spec2max = [(i+j) for i,j in zip(spec2,uncert2)]
spec2min = [(i-j) for i,j in zip(spec2,uncert2)]

# path to all bayesian unfolding data 
files_json = 'files.json'

#number of parameters and their names and init guesses for bayes
nparam = 4
param_names = [#r'$\mu_{fast,1}$',
               r'$\sigma_{fast,1}$',r'$C_{fast,1}$',
               #r'$\sigma_{fast,2}$',r'$C_{fast,2}$',
               #r'$\sigma_{fast,3}$',r'$C_{fast,3}$',
               r'$T_W$', r'$C_W$'
               ]
guesses = [#15.5-2.0945, 
           0.35, 5e5, 
           #0.5, 1e5,
           #0.5, 1e4,
           1,  1e6
           ]

# initialise object to load bayes model
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

# find chi squared for both spectra
chisq_spec1 = unfold.get_chi_squared(spec1)
chisq_spec2 = unfold.get_chi_squared(spec2)
print('Chi squared =',chisq_spec1, 'and', chisq_spec2)

# find mean spectrum uncertainty
std_spec1 = 100*np.mean(uncert1/spec1)
std_spec2 = 100*np.mean(uncert2/spec2)
print('mean spectrum std =',std_spec1, 'and', std_spec2, 'pc')

# group structure shorthand
gs = unfold.group_structure

# plot all solution spectra
fig, ax = plt.subplots(1,figsize=(10,6))
ax.step(gs, spec1, label='Bayesian',c='m',where='pre')
ax.fill_between(gs, spec1min,spec1max,
                alpha=0.2,step='pre',color='m')
ax.step(gs, spec2,label='GRAVEL',c='g',where='pre')
ax.fill_between(gs, spec2min,spec2max,
                alpha=0.2,step='pre',color='g')
#ax.set_xscale('log')
ax.set_xlim(0,20)
ax.set_xlabel('Neutron energy (MeV)')
ax.set_yscale('log')
ax.set_ylim(1e1)
ax.set_ylabel('Flux per energy bin (n cm$^{-2}$ s$^{-1}$)')
ax.grid()
ax.legend(loc="lower left",frameon=True, fontsize=18,
          fancybox=False,facecolor='white')
fig.tight_layout()
plt.savefig(f'combo_spectrum2.png')
