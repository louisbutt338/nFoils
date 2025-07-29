from nfoils.reaction import ResponseFunction
import sandy
import os

# get NJOY response functions - need internet for sandy to get library data 

# set NJOY env variable if you haven't already
os.environ['NJOY'] = '/Users/ljb841@student.bham.ac.uk/NJOY2016/bin/njoy'

# set energy grids (feel free to make own in sandys format)
ek=sandy.energy_grids.SCALE238

# set library (check if availble in sandy)
library = 'jeff_33' # endfb_71 endfb_80 jendl_40u jeff_33  tendl_21

# filename for input data
data_file_name = 'foil_data'

# list of reaction labels for plotting (must match data in input file)
reaction_labels = [r'${}^{115}$In(n,$\gamma$)',
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

# run that shi
response_fn_plot = ResponseFunction(ek,library,data_file_name,reaction_labels)
response_fn_plot.run_xs()