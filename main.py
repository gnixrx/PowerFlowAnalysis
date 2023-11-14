import numpy as np
from data_reader import DataReader
from data_writer import DataWriter
from nodal_analysis import NodalAnalysis
from power_analysis import PowerAnalysis

def main():
    # Set up base case
    data_read = DataReader("system_basecase.xlsx")
    nodal = NodalAnalysis(data_read.line_data)
    power = PowerAnalysis(data_read.bus_data, nodal.get_y_matrix())
    if (nodal.get_node_max() != power.get_node_max()):
        print("Error: Maximum nodes in line data does not match bus data.")
        exit()

    data_write = DataWriter()
    data_write.write_y_matrix(nodal.get_y_matrix())

    # Run base case.
    print("Running base case.")
    print(power.update())

    print(data_read.bus_data)
    print()

    # Set up contingency case 1
    # One circuit of line 1 <-> 2 are taken out of service.
    # data.line_data[0].set_state("Half")
    # print(data.line_data[0])

    # Set up contingency case 2
    # Lines 4 <-> 5, 8 <-> 9, 10 <-> 11 are taken out of service.
    # data.line_data[0].set_state("On")
    # for n in [6, 14, 15]:
    #     data.line_data[n].set_state("Off")
    #     print(data.line_data[n])

if __name__ == "__main__":
    main()