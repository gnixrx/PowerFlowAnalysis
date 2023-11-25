import numpy as np
from global_setting import s_base
from data_bus import BusData
from scipy.sparse import csr_array

"""
The LineData class describes the connections between nodes as physical lines in the pi model, or as a 
    simple representation of a transformer.
:param going_from: The node identifer where this connection is being made from.
:type going_from: int
:param going_to: The node identifier where this connection is being made to.
:type going_to: int
:param resistance: The series resistance within the pi model of a transmission line.
:type resistance: float
:param reactance: The series reactance within the pi model of the transmission line or the reactance in 
    the simple representation of a transformer.
:type reactance: float
:param total_shunt_susceptance: The total amount of susceptance in the pi model of the transmission line.
:type total_shunt_susceptance: float
:param maximum_capacity: Limit of power transmission capacity in MVA.
:type maximum_capacity: int
"""
class LineData:
    ### ----------------------------------------- Getters -----------------------------------------
    @property
    def bus_from(self):
        return self._f

    @property
    def bus_to(self):
        return self._t

    @property
    def p_from(self):
        p = self.__line_flow_calc(self._f, self._t)
        if p != "Unknown":
            p = np.real(p)
        return p

    @property
    def q_from(self):
        q = self.__line_flow_calc(self._f, self._t)
        if q != "Unknown":
            q = np.imag(q)
        return q

    @property
    def p_to(self):
        p = self.__line_flow_calc(self._t, self._f)
        if p != "Unknown":
            p = np.real(p)
        return p

    @property
    def q_to(self):
        q = self.__line_flow_calc(self._t, self._f)
        if q != "Unknown":
            q = np.imag(q)
        return q

    @property
    def overloaded(self):
        """
        Returns the overload status of the line.
        :return: Overload status "Unknown", True, or False
        :rtype: string or boolean
        """
        # Don't bother if power is unknown.
        if (self.p_from == "Unknown" or self.q_from == "Unknown"
                or self.p_to == "Unknown" or self.q_to == "Unknown"):
            return "Unknown"

        # Find power
        from_mva = np.sqrt(self.p_from ** 2 + self.q_from ** 2)
        to_mva = np.sqrt(self.p_to ** 2 + self.q_to ** 2)

        # Set multipler for state.
        if self._state == "Off": m = 0
        elif self._state == "Half":  m = 0.5
        else: m = 1

        if from_mva <= m * self._fmax or to_mva <= m * self._fmax:
            return False
        return True


    @property
    def state(self):
        """
        Return the current state of this bus.
        :return: Bus state
        :rtype: string
        """
        return self._state


    @property
    def type(self):
        """
        This function returns the type of connection between the nodes as either a transformer,
            or a line depending on if the connection has a shunt susceptance.
        :return: Type of connection being made between the nodes.
        :rtype: string
        """
        return "Transformer" if self._half_b == 0 else "Line"


    @property
    def max_node(self):
        """
        Returns the maximum node referenced by this line either its to or from attribute
        :return: Maximum identifier of the to or from attribute
        :rtype: int
        """
        return max(self._f.id, self._t.id)

    ### ------------------------------------------ Setters ------------------------------------------

    @state.setter
    def state(self, new_state: str):
        """
        Sets up the state of the line.
        :param new_state: New state of "Off" or "Half" otherwise will default to the on state.
        :type new_state: string
        :return: None
        :rtype: None
        """
        if new_state == "Off":
            self._state = "Off"
        elif new_state == "Half":
            self._state = "Half"
        else:
            self._state = "On"

    ### ----------------------------------------- Functions -----------------------------------------

    def y_stamp(self, size):
        """
        Returns the linedata stamp for the Y (admittance) matrix
        :param size: The shape size of the NodalAnalysis Y matrix
        :type size: int
        :return: The sparse array for the conductance in this connection between nodes in the network.
        :rtype: class: csr_array
        """
        # Set multipler for state.
        if self._state == "Off": m = 0
        elif self._state == "Half": m = 0.5
        else: m = 1

        # Create the stamp values for diagonal and nondiagonals for this line.
        y_diag_stamp = m * (1 / (self._r + 1j * self._x) + 1j * self._half_b)
        y_nondiag_stamp = m * -1 / (self._r + 1j * self._x)

        # Create and return the sparse array.
        row = [self._f.id, self._t.id, self._f.id, self._t.id]
        col = [self._f.id, self._t.id, self._t.id, self._f.id]
        data = [y_diag_stamp, y_diag_stamp, y_nondiag_stamp, y_nondiag_stamp]
        return csr_array((data, (row, col)), shape=(size, size))


    def __line_flow_calc(self, k, i):
        """
        Private function to calculate the line flow from bus to bus. Only up to date after
            power_analysis.update()
        :param k: From bus
        :type k: BusData
        :param i: To bus
        :type i: BusData
        :return: Complex power
        :rtype: complex float
        """
        # Unknown Voltages at each bus
        if k.V == None or i.V == None:
            return "Unknown"

        # Set multipler for state. Y_total = Y_1 + Y_2, if Y_1 = Y_2 = Y, then Y = Y_total * 0.5
        if self._state == "Off": m = 0
        elif self._state == "Half": m = 0.5
        else: m = 1

        Vk = k.V * np.cos(k.Th) + 1j * k.V * np.sin(k.Th)
        Vi = i.V * np.cos(i.Th) + 1j * i.V * np.sin(i.Th)

        I_ki = m * ((Vk - Vi) / (self._r + 1j * self._x) + 1j * self._half_b * Vk)
        result = Vk * np.conjugate(I_ki)
        return result


    def __init__(self,
                 going_from: BusData,
                 going_to: BusData,
                 resistance: float,
                 reactance: float,
                 total_shunt_susceptance: float,
                 maximum_capacity: int):
        self._f = going_from # Connection going from BusData
        self._t = going_to # Connection going to BusData
        self._r = resistance # Resistance of connection in pu
        self._x = reactance # Reactance of connection in pu
        self._half_b = total_shunt_susceptance / 2 # Half the total shunt susceptance in pu
        self._fmax = maximum_capacity / s_base # Maximum transmission capacity in pu
        self._state = "On" # Sets state to on for this line


    def __repr__(self):
        return (f"{self.__class__.__name__}> Connection: {self._f.id + 1} <-> {self._t.id + 1}, "
                f"Type: {self.type}, State: {self.state}")