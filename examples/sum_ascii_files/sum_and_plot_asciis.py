""" example for summing ascii files together and plotting the result
"""
from nfoils.ascii import AsciiSummer

# path to dir with asciis in
folder_path = 'example_asciis'

# name of files without their numbers at the end 
ascii_filetag = 'uBB_20s_x60_20mins'

# list of file numbers to sum together i.e. [0,3] = '000' to '003'
file_numbers = [0,3]

# sum the asciis and plot
ascii_summer = AsciiSummer(folder_path,ascii_filetag,file_numbers)
ascii_summer.run()