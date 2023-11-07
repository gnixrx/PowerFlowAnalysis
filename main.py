import numpy as np
from data_reader import DataReader
from nodal_analysis import NodalAnalysis
from global_setting import s_base

def main():
    np.set_printoptions(suppress=True)
    data = DataReader()

    nodal = NodalAnalysis(data.line_data, data.bus_data)
    print(np.round(nodal.get_y_matrix(), 2))


if __name__ == "__main__":
    main()