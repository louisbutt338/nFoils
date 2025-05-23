from nfoils.fispact_tools import JsonPlotter

# the main worker
def main():
    # input directory with the json file in (which is also the directory where the plots will be saved)
    directory = '/Users/ljb841@student.bham.ac.uk/fispact/WORKSHOP/demo_hcll/1hr1g_fullflux_fe'
    # input json filename (not extension)
    json_name = 'hcll_fw_fe'

    json_plot = JsonPlotter(directory,json_name) 
    json_plot.run()

if __name__ == '__main__':
    main() 


