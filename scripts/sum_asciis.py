from nfoils.ascii_summer import AsciiSummer

# the main worker
def main():
    ascii_filetag = 'ubb_400s_x60_6dot7hrs'
    first_file_number = 0
    last_file_number = 55
    folder_path = '/Users/ljb841@student.bham.ac.uk/gamma_spec/deuteron_hpge/hpge_results_g11_291124/long-lived/y'
    ascii_summer = AsciiSummer(folder_path,ascii_filetag,first_file_number,last_file_number)
    ascii_summer.write_ascii()

if __name__ == '__main__':
    main() 


