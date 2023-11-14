import numpy as np
import pandas as pd
from data_line import LineData
from data_bus import BusData

class DataWriter():
    def write_y_matrix(self, y_mat):
        # Create data frame
        l = np.shape(y_mat)[0] + 1
        c_i = np.char.mod('%d', np.arange(1, l))
        g = np.real(y_mat.toarray())
        b = np.imag(y_mat.toarray())
        matrix = y_mat.toarray()
        self._y_mat_df = pd.DataFrame(matrix, columns=c_i, index=c_i)

        # Write sheet
        with self._xls as writer:
            self._y_mat_df.to_excel(writer, sheet_name="Acceptance Matrix", float_format="%.3f")

    # Initialize the data writer class.
    def __init__(self, file="system_output.xlsx"):
        # Initialize writer
        self._xls = pd.ExcelWriter(file)

