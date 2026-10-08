""" example for postprocessing spectrum after unfolding
"""
from nfoils.unfold import BayesianUnfolding
import numpy as np
import json

# linear group structure number
gs_num = 51

# path to files json containing paths for all the data
files_json = 'files.json'

# edit the files json file in this script if needed
with open(files_json, 'r') as f:
    data = json.load(f)
    data["response_matrix"]= f"data_linear_gs/for_unfolding_{gs_num}/response_matrix.txt"
    data["response_matrix_uncerts"]= f"data_linear_gs/for_unfolding_{gs_num}/response_matrix_uncertainties.txt"
    data["group_structure"] = f"data_linear_gs/for_unfolding_{gs_num}/group_structure.txt"
with open('files.json','w') as f:
    json.dump(data,f,indent=4)

# path to results folder
#results_path = f'good_lhood_results/diff_gs_5param_weiss_260226/lin{gs_num}'
results_path = 'good_lhood_results/diff_priors_lin51_240226/1freepeak_maxwell'

# how many parameters in the model
nparam = 5

# names of parameters
param_names = [r'$\mu_{fast,1}$',
               r'$\sigma_{fast,1}$',r'$C_{fast,1}$',
               #r'$\sigma_{fast,2}$',r'$C_{fast,2}$',
               #r'$\sigma_{fast,3}$',r'$C_{fast,3}$',
               r'$T_W$', r'$C_W$'
               ]

# set initial guesses for each parameter
guesses = [15.5-2.0945, 
           0.35, 5e5, 
           #0.5, 1e5,
           #0.5, 1e4,
           1,  1e6
           ]

# initialise object to load data
unfold = BayesianUnfolding(files_json,nparam,param_names,guesses)

# create neutron flux model for plotting the parameters
def model(theta,energy):

    # unpack the tuple of parameters
    (mean1,
     sigma1,peak1,
     #sigma2,peak2,
     #sigma3,peak3,
     mode4,peak4,
     #alpha5,scale5
     )=theta

    # create the model
    #if energy>1:
    flux = (unfold.gaussian(mean1,sigma1,peak1,energy)+
            #unfold.gaussian(15.5-6.45,sigma2,peak2,energy)+
            #unfold.gaussian(15.5-8.85,sigma3,peak3,energy)+
            unfold.maxwellian(mode4,peak4,energy)
            )
    #else:
    #    flux = (unfold.weisskopf(mode4,peak4,energy)+
    #            unfold.powerlaw(alpha5,scale5,energy))
    return flux

# load calculated parameters and standard deviations from results file
#params_fast = np.loadtxt('results_test.txt', usecols=(1,))
#params_slow = np.loadtxt('results_slow.txt', usecols=(1,))
#params = np.append(params_fast,params_slow)
params = np.loadtxt(f'{results_path}/params.txt')

# get covariance matrix
#covs_fast = np.loadtxt('results_test.txt', usecols=(2,))
#covs_slow = np.loadtxt('results_slow.txt', usecols=(2,))
#covs = np.append(stds_fast,stds_slow)
cov_matrix = np.loadtxt(f'{results_path}/covs.txt')

# sample cov matrix by cholesky decomp
chol = np.linalg.cholesky(cov_matrix)
rand_samples = np.random.randn(100000000,nparam)
true_samples = rand_samples @ chol.T
true_cov_matrix = np.cov(true_samples.T)

# calc stds from new cov matrix
stds=np.diag(np.sqrt(true_cov_matrix))

# unphysical values to cutoff from end of spectrum
cutoff = None

# find the mean param uncertainty
print(f'mean param std = {100*np.mean(stds/params)} pc')

# get the different spectra out from the params and cov matrix
all_spectra = unfold.get_spectra(model,params,true_cov_matrix,cutoff)

# find the mean spectrum uncertainty
mean_std_spec = 100*np.mean((all_spectra[1]-all_spectra[0])/all_spectra[0])
print(f'mean spectrum std = {mean_std_spec} pc')

# plot and dump spectrum
#unfold.plot_simple_spectrum(*all_spectra,f'{results_path}/spectrum4')
#unfold.dump_spectrum(all_spectra[0],all_spectra[1],f'{results_path}/spectrum4')

#and reduced chi squared
chi_squared = unfold.get_chi_squared(all_spectra[0])
print('Chi squared =',chi_squared)