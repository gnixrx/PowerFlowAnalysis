from global_setting import s_base
from scipy.sparse import csr_array

"""
The LineData class describes the connections between nodes as physical lines in the pi model, or as a simple 
    representation of a transformer.
:param going_from: The node identifer where this connection is being made from.
:type going_from: int
:param going_to: The node identifier where this connection is being made to.
:type going_to: int
:param resistance: The series resistance within the pi model of a transmission line.
:type resistance: float
:param reactance: The series reactance within the pi model of the transmission line or the reactance in the 
    simple representation of a transformer.
:type reactance: float
:param total_shunt_susceptance: The total amount of susceptance in the pi model of the transmission line.
:type total_shunt_susceptance: float
:param maximum_capacity: Limit of power transmission capacity in MVA.
:type maximum_capacity: int
"""
class LineData:
    def set_state(self, new_state: str):
        if new_state == "Off":
            self._state = "Off"
        elif new_state == "Half":
            self._state = "Half"
        else:
            status = "On"

    def type(self):
        """
        This function returns the type of connection between the nodes as either a transformer, or a line depending on
        if the connection has a shunt susceptance.
        :return: Type of connection being made between the nodes.
        :rtype: string
        """
        return "Transformer" if self._half_b == 0 else "Line"

    def max_node(self):
        """
        Returns the maximum node referenced by this line either its to or from attribute
        :return: Maximum identifier of the to or from attribute
        :rtype: int
        """
        return max(self._f, self._t)

    # Return linedata stamps in pi model (where transformers are only a reactance)

    def y_stamp(self, size):
        """
        Returns the linedata stamp for the Y (admittance) matrix
        :param size: The shape size of the NodalAnalysis Y matrix
        :type size: int
        :return: The sparse array for the conductance in this connection between nodes in the network.
        :rtype: class: csr_array
        """
        # Set multipler for state
        if self._state == "Off": m = 0
        elif self._state == "Half": m = 0.5 # Y_total = Y_1 + Y_2 therefore Y_2 = Y_total - Y_1
        else: m = 1

        y_diag_stamp = m * 1 / (self._r + 1j * self._x) + 1j * self._half_b
        y_nondiag_stamp = m * -1 / (self._r + 1j * self._x)

        row = [self._f, self._t, self._f, self._t]
        col = [self._f, self._t, self._t, self._f]
        data = [y_diag_stamp, y_diag_stamp, y_nondiag_stamp, y_nondiag_stamp]
        return csr_array((data, (row, col)), shape=(size, size))

    # Initialize the class
    def __init__(self, going_from: int, going_to: int, resistance: float, reactance: float, total_shunt_susceptance: float, maximum_capacity: int):
        self._f = going_from - 1 # connection going from index
        self._t = going_to - 1 # connection going to index
        self._r = resistance # resistance of connection in pu
        self._x = reactance # reactance of connection in pu
        self._half_b = total_shunt_susceptance / 2 # half the total shunt susceptance in pu
        self._fmax = maximum_capacity / s_base # maximum transmission capacity in pu
        self._state = "On"

    def __repr__(self):
        return f"{self.__class__.__name__}> Connection: {self._f + 1} <-> {self._t + 1}, Type: {self.type()}, State: {self._state}"