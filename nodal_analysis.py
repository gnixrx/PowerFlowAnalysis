import numpy as np
from scipy.sparse import csr_array
from data_line import LineData
from data_bus import BusData

"""
NodalAnalysis class describes the circuit network of lines and connections of busses of the power system.
:param line_data: A list of LineData which represent the connections between nodes in the graph of the network in the power system. Such as
                  lines or transformers connecting the different busses in the power system.
:type line_data: list of class: LineData
:param bus_data: A list of BusData which represent the attributes of the nodes in the graph of the network in the power system, including
                 power consumed and generated at the node.
:type bus_data: list of class: BusData

"""
class NodalAnalysis:
    def get_gb_matrix(self):
        """
        Returns the separated conductance and susceptance matrices of the current network.
        :return: The separataed _Y matrix as the conductance matrix _G, and the susceptane matrix _B separately.
        :rtype: array, array
        """
        Y = self.get_y_matrix()
        return np.real(Y), np.imag(Y)

    # Return the combined Y admittance matrix as array
    def get_y_matrix(self):
        """
        Returns the combined complex admittance matrix of the current network.
        :return: The combined _Y matrix
        :rtype: array
        """

        # Create admittance matrix
        self._Y = csr_array((self._node_max, self._node_max), dtype=float) # Y Admittance Matrix
        for line in self._line_data:
            self._Y = self._Y + line.y_stamp(self._node_max)

        return self._Y.toarray()

    # Initialize the class
    def __init__(self, line_data: LineData, bus_data: BusData):
        self._line_data = line_data

        # Create node count
        ld_max = max(line.max_node() for line in line_data)
        bd_max = bus_data.size
        self._node_max = max(ld_max, bd_max) # Maximum number of nodes in the system


    def __repr__(self):
        return f"{self.__class__.__name__}> Number of nodes: {self.nodes_max}"