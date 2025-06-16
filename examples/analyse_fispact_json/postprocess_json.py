from nfoils.fispact_tools import JsonPlotter

# the main worker
def main():
    # input working directory with the json file in and to save plot
    directory = ''
    # input json filename (not extension)
    json_name = 'hcll_fw_fe'

    json_plot = JsonPlotter(directory,json_name) 
    json_plot.run()

if __name__ == '__main__':
    main() 


