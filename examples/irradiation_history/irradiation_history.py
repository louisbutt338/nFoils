""" example for converting time bin data into fispact irrad history
"""

from nfoils.history import IrradTimings

# path to the rate data txt file with
# time(s) and count rate (1/s) in the two columns
filename = "rate"

# desired name of the txt output file
output_filename = "irrad_history"

# init
history = IrradTimings(filename, output_filename)

# enter conversion factor for going from rate data cps 
# to number of charged particles per sec
# here, example based on 12uA=190k counts
cps_to_particles_per_s = 6.24151e12 * (12/190000)
# then write a fispact history file
history.fispact_hist_writer(cps_to_particles_per_s)

# OR enter conversion factor for rate to current
# again example based on 12uA=190k counts
rate_to_current = (12/190000)
# print an array of currents 
history.current_printer(rate_to_current)
