""" example for summing ascii files together and plotting the result
"""
from nfoils.ascii import AsciiSummer

# path to dir with asciis in
folder_path = ('/Users/ljb841@student.bham.ac.uk/gamma_spec/data/dli_sep25/'
               'short-lived/al_1cm')

# name of files without the '_xxx' numbers at the end 
ascii_filetag = 'louis_dLi_20s_x60_20mins'

# list of file numbers to sum together i.e. [0,3] = '000' to '003'
file_numbers = [0,59]

# sum the asciis and plot
ascii_summer = AsciiSummer(folder_path,ascii_filetag,file_numbers)
ascii_summer.run()