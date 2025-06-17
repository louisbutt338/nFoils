from nfoils.ascii_summer import AsciiSummer

folder_path = 'fe_spectra'
#name of files withoiut their number at the end 
ascii_filetag = 'fe_ubb_280324'
first_file_number = 0
last_file_number = 0

ascii_summer = AsciiSummer(folder_path,ascii_filetag,first_file_number,last_file_number)
ascii_summer.run()