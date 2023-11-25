import string
import numpy as np
import pandas as pd
from global_setting import s_base
from data_line import LineData
from data_bus import BusData

"""
Data writer class creates the excel outputs for the analysis.
:param file: The path string of the output file. Defaults to system_output.xlsx.
:type file: string
"""
class DataWriter():
    ### ------------------------------------------ Getters ------------------------------------------
    @property
    def file(self):
        """
        Returns the file name.
        :return: File name.
        :rtype: string
        """
        return self._file

    ### ----------------------------------------- Functions -----------------------------------------
    def add_y_matrix(self, y_mat):
        """
        Adds an acceptance (Y) matrix to the write queue.
        :param y_mat: Sparse array of the acceptance matrix.
        :type y_mat: csr_array
        :return: None
        :rtype: none
        """
        # Create data frame
        l = np.shape(y_mat)[0] + 1
        c_i = np.char.mod('%d', np.arange(1, l))
        y = np.round(y_mat.toarray(), 3)
        self._y_mat_df.append(pd.DataFrame(y, columns=c_i, index=c_i))

    def add_power_iterations(self, iterations, exec_time, mm_record):
        """
        Adds power iterations data to the write queue.
        :param iterations: Number of iterations the power calculation went through.
        :type iterations: int
        :param exec_time: Time it took in seconds to calculate power iterations.
        :type exec_time: float
        :param mm_record: Record of the maximum mismatches for each iteration for both active and
            reactive power
        :type mm_record: list
        :return: None
        :rtype: none
        """
        col = ['P Bus', 'P Mismatch', 'Q Bus', 'Q Mismatch']
        for row in mm_record:
            row[0] = row[0].id + 1
            row[2] = row[2].id + 1
        self._mm_record_df.append(pd.DataFrame(mm_record, columns=col))

    def add_bus_line_result(self, bus_data, line_data):
        """
        Adds the bus and line results to the outputs
        :param bus_data: The list of bus data classes which have been updated by the
            poweranalysis.update function.
        :type bus_data: list of BusData
        :param line_data: The list of line data classes
        :type line_data: list of LineData
        :return: None
        :rtype: none
        """
        # Create bus results
        col = ['Bus', 'Type', 'V (p.u.)', 'Angle (degrees)', 'Real Power Gen (MW)',
               'Reactive Power Gen (MVAr)', 'Real Power Load (MW)', 'Reactive Power Load (MVAr)',
               'Within Voltage Limits' ]
        b_result = []
        for bus in bus_data:
            b_result.append([bus.id + 1, bus.type, bus.V, np.degrees(bus.Th),
                             bus.P_gen * s_base, bus.Q_gen * s_base,
                             bus.P_load * s_base, bus.Q_load * s_base,
                             str(bus.in_V_limit)])

        self._bus_result_df.append(pd.DataFrame(b_result, columns=col))

        # Create line results
        col = ['Bus k', 'Bus i', 'Real Power Flow k to i (MW)', 'Reactive Power Flow k to i (MVAr)',
               'Real Power Flow i to k (MW)', 'Reactive Power Flow i to k (MVar)',
               'Real Power Losses (MW)', 'Reactive Power Losses (MVar)', 'Overload']

        l_result = []
        for line in line_data:
            l_result.append([line.bus_from.id + 1, line.bus_to.id + 1,
                             line.p_from * s_base, line.q_from * s_base,
                             line.p_to * s_base, line.q_to * s_base,
                             abs(line.p_from + line.p_to) * s_base,
                             abs(line.q_from + line.q_to) * s_base,
                             str(line.overloaded)])

        self._line_result_df.append(pd.DataFrame(l_result, columns=col))

    def write(self):
        """
        Writes out the excel output file.
        :return: Success
        :rtype: boolean
        """
        # Write sheet
        with self._xls as writer:
            # Write Y Matrix
            for i in range(len(self._y_mat_df)):
                s_name = "Y Matrix"
                if len(self._y_mat_df) > 1:
                    s_name = f"{s_name} {i + 1}"
                self._y_mat_df[i].to_excel(writer, sheet_name=s_name)
            # Write Convergence
            for i in range(len(self._mm_record_df)):
                s_name = "Convergence Record"
                if len(self._mm_record_df) > 1:
                    s_name = f"{s_name} {i + 1}"
                self._mm_record_df[i].to_excel(writer, sheet_name=s_name)
            # Write Bus/Line Results
            for i in range(len(self._bus_result_df)):
                s_name = "Bus Result"
                if len(self._bus_result_df) > 1:
                    s_name = f"{s_name} {i + 1}"
                self._bus_result_df[i].to_excel(writer, sheet_name=s_name, index=False)

                s_name = "Line Result"
                if len(self._line_result_df) > 1:
                    s_name = f"{s_name} {i + 1}"
                self._line_result_df[i].to_excel(writer, sheet_name=s_name)

            return True
        return False


    # Initialize the data writer class.
    def __init__(self, file="system_output.xlsx"):
        self._file = file
        self._xls = pd.ExcelWriter(file, engine='openpyxl')
        self._y_mat_df = []
        self._mm_record_df = []
        self._bus_result_df = []
        self._line_result_df = []

