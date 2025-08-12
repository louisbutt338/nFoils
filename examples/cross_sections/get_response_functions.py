from nfoils.reaction import PostprocessReactions
from nfoils.reaction import IsotopicSpectrumUncertainty
import sandy
import os

# set NJOY env variable if you haven't already
os.environ['NJOY'] = '/Users/ljb841@student.bham.ac.uk/NJOY2016/bin/njoy'

# set energy grids (feel free to make own in sandys format)
ek=sandy.energy_grids.SCALE238

# set library (check if availble in sandy) - need internet to get library data 
library = 'jeff_33' # endfb_71 endfb_80 jendl_40u jeff_33  tendl_21

# filename for input data
data_file_name = 'foil_data'

# list of reaction labels for plotting (must match data in input file)
reaction_labels = [r"${}^{115}$In(n,n')", 
                   r'${}^{65}$Cu(n,p) *',
                   r'${}^{56}$Fe(n,p)']

# get some response functions and nuclear data uncertainties
reactions = PostprocessReactions(ek,library)
reactions.run_rf(data_file_name,reaction_labels)
#reactions.run_stdev(data_file_name,reaction_labels)

# or get some reaction-dependent spectrum uncertainties 
#uncertainties = IsotopicSpectrumUncertainty(ek,library)
#spectrum_file = 'spectra'
#uncertainties.get_isotopic_uncertainties(spectrum_file,data_file_name)