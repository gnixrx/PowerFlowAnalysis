import numpy as np
from scipy.sparse import csr_array
from data_line import LineData

"""
NodalAnalysis class describes the circuit network of lines and connections of busses of the power system.
:param line_data: A list of LineData which represent the connections between nodes in the graph of the network in the power system. Such as
                  lines or transformers connecting the different busses in the power system.
:type line_data: np.Array of LineData

"""
class NodalAnalysis:
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

        return self._Y

    def get_node_max(self):
        """
        Returns the maximum node numbers in the network from the line data
        :return: Maximum node
        :rtype: int
        """
        return self._node_max

    # Initialize the class
    def __init__(self, line_data: np.array):
        self._line_data = line_data

        # Create node count
        self._node_max = max([line.max_node() for line in line_data]) + 1 # Maximum node refered in each line


    def __repr__(self):
        return f"{self.__class__.__name__}> Number of nodes: {self._node_max}"