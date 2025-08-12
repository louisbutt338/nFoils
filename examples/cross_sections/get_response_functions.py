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

#reactions = PostprocessReactions(ek,library,data_file_name,reaction_labels)
# get some response functions 
#reactions.run_rf()
# or some uncertainties
#reactions.run_stdev()

# get some specific spectrum uncertainties 
spectrum_file = 'spectra'
uncertainties = IsotopicSpectrumUncertainty(ek,library,data_file_name,reaction_labels)
uncertainties.get_isotopic_uncertainties(spectrum_file)