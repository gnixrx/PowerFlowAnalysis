import numpy as np
from data_reader import DataReader
from data_writer import DataWriter
from nodal_analysis import NodalAnalysis
from power_analysis import PowerAnalysis

def main():
    # Initialize
    data_read = DataReader("system_basecase.xlsx")
    data_write = DataWriter()
    nodal = NodalAnalysis(data_read.line_data)

    data_write.add_y_matrix(nodal.get_y_matrix())  # Deliverable #1
    power = PowerAnalysis(data_read.bus_data, nodal.get_y_matrix())

    # Run base case.
    print("Running base case.")
    print(power.update())
    data_write.add_power_iterations(*power.update_data) # Deliverable #2
    data_write.add_bus_line_result(data_read.bus_data, data_read.line_data) # Deliverable #3

    # Set up contingency case 1
    # One circuit of line 1 <-> 2 are taken out of service.
    # print("Running contingency case 1.")
    # data_read.line_data[0].state = "Half"
    # print(data_read.line_data[0])
    # print(power.update())
    # data_write.add_bus_line_result(data_read.bus_data, data_read.line_data) # Deliverable #4

    # Set up contingency case 2
    # Lines 4 <-> 5, 8 <-> 9, 10 <-> 11 are taken out of service.
    # print("Running contingency case 2.")
    # data.line_data[0].state = "On"
    # for n in [6, 14, 15]:
    #     data.line_data[n].state = "Off"
    #     print(data.line_data[n])
    # print(power.update())
    # data_write.add_bus_line_result(data_read.bus_data, data_read.line_data) # Deliverable #5

    # Create excel deliverable
    data_write.write()

    print()

if __name__ == "__main__":
    main()