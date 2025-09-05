""" example for converting time bin data into fispact irrad history
"""

from nfoils.history import IrradTimings

filename = "rate.txt"
output_filename = "irrad_history.txt"

# acquisition start time was 29/11/24 13:46:35
# txt data is time(s) and count rate (1/s)

#run 
fispact_history = IrradTimings(filename, output_filename)
fispact_history.run()
