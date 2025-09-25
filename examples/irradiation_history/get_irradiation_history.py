""" example for converting time bin data into a fispact neutron irrad history 
and an array of cyclotron charged particle currents using conversion factors
"""

from nfoils.history import IrradTimings

# path to the rate data txt file with
# time(s) and count rate (1/s) in the two columns
filename = "rate"

# desired name of the txt output file
output_filename = "irrad_history"

# initialise
history = IrradTimings(filename, output_filename)

# enter conversion factors for going from:
# rate data cps --> charged particle current uA --> charged particles/s
# here, example based on 12uA=190k counts
cps_to_current = (12/190000)
cps_to_particles_per_s = 6.24151e12 * cps_to_current

# print an array of currents 
history.current_printer(cps_to_current)

# write a fispact history file
history.fispact_hist_writer(cps_to_particles_per_s)
