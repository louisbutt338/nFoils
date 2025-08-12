from nfoils.ascii import AsciiSummer

folder_path = 'example_asciis'
#name of files withoiut their number at the end 
ascii_filetag = 'uBB_20s_x60_20mins'
first_file_number = 0
last_file_number = 3

ascii_summer = AsciiSummer(folder_path,ascii_filetag,first_file_number,last_file_number)
ascii_summer.run()