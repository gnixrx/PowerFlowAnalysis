import numpy as np
from data_reader import DataReader
from nodal_analysis import NodalAnalysis
from power_analysis import PowerAnalysis

def main():
    np.set_printoptions(suppress=True)

    # Set up base case
    data = DataReader("simple_basecase.xlsx")
    nodal = NodalAnalysis(data.line_data)
    power = PowerAnalysis(data.bus_data, nodal.get_y_matrix())
    if (nodal.get_node_max() != power.get_node_max()):
        print("Error: Maximum nodes in line data does not match bus data.")
        exit()

    power.update()

    print(data.bus_data)

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