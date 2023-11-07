import pandas as pd
import numpy as np
from data_line import LineData
from data_bus import BusData

"""
Loads and organizes the provided specification for the power system to be analysed from an excel file in a similar format 
to the system_basecase.xlsx file provided.
:param file: The filename describing the power system to be analysed.
:type file: string
"""
class DataReader:
    # Getters
    @property
    def line_data(self):
        """
        Returns a list of the LineData types.
        :return: Returns a simple python list of the LineData described in the specification excel file.
        :rtype: list of LineData
        """
        return self._line_data
    @property
    def bus_data(self):
        """
        Returns a list of the BusData types.
        :return: Returns a simple python list of the BusData described in the specification excel file.
        :rtype: list of BusData
        """
        return self._bus_data

    # Initialize the data class.
    def __init__(self, file="system_basecase.xlsx"):
        # Check to see if the system excel definition is present
        if file[-5:] == ".xlsx" or file[-4:] == ".xls":
            self.file = file
        else:
            print("System definition file is not present or an excel file. Reverting to default.")

        line_data_raw = pd.read_excel(open(self.file, "rb"), sheet_name="LineData")
        bus_data_raw = pd.read_excel(open(self.file, "rb"), sheet_name="BusData")

        # Check to see if the data is populated
        if line_data_raw.size == 0:
            print("Line data not present. Reconcider system definition file.")

        if bus_data_raw.size == 0:
            print("Bus data no present. Reconcider system definition file.")

        # Populate the data as objects
        # List of LineData
        self._line_data = np.array([ LineData(row['From'], row['To'], row['Rtotal, p.u.'], row['Xtotal, p.u.'], row['Btotal, p.u.'], row['Fmax, MVA'])
                                     for row in line_data_raw.to_dict(orient='records') ])
        # List of BusData
        self._bus_data = np.array([ BusData(row['Bus #'], row['P MW'], row['Q MVAr'], row['Type'], row['P Gen'], row['V Set'])
                                    for row in bus_data_raw.to_dict(orient='records') ])

    def __repr__(self):
        return f'{self.__class__.__name__}> File: {self.file}'
