import numpy as np
from global_setting import s_base

"""
The BusData class describes the attributes of the busses in the power system.
:param b: The node/bus identifer.
:type going_from: int
:param p: Load active power loss in MW.
:type resistance: float
:param q: Load reactive power loss in MVAr.
:type reactance: float
:param type: Bus type 'S' for Slack bus, 'G' for PV bus, 'D' for PQ bus.
:type type: chr
:param p_gen: Generated active power in MW.
:type p_gen: int
:param v_set: Voltage set in pu.
:type v_set: float
"""

class BusData:
    ### ----------------------------------------- Getters -----------------------------------------
    @property
    def id(self):
        """
        Returns the bus number
        :return: Bus Number
        :rtype: int
        """
        return self._b

    @property
    def P(self):
        """
        Returns the active power at this bus.
        :return: Active power
        :rtype: float
        """
        return self._P

    @property
    def P_gen(self):
        """
        Returns how much power is generated at this bus.
        :return: Active power generated at this bus.
        :rtype: float
        """
        return self.P + self._p_load

    @property
    def P_load(self):
        """
        Returns the load at this bus.
        :return: Active power load at this bus.
        :rtype: float
        """
        return self._p_load

    @property
    def Q(self):
        """
        Returns the reactive power at this bus.
        :return: Reactive power
        :rtype: float
        """
        return self._Q

    @property
    def Q_gen(self):
        """
        Returns the reactive power generated at this bus.
        :return: Reactive power
        :rtype: float
        """
        return self.Q + self._q_load

    @property
    def Q_load(self):
        """
        Returns the reactive power load at this bus.
        :return: Reactive power
        :rtype: float
        """
        return self._q_load

    @property
    def V(self):
        """
        Returns the voltage at this bus in pu.
        :return: Voltage
        :rtype: float
        """
        return self._V

    @property
    def Th(self):
        """
        Returns the voltage angle at this bus in radians.
        :return: Voltage angle, theta
        :rtype: float
        """
        return self._Th

    @property
    def type(self):
        """
        Returns the type of this bus, Slack, PV, or PQ.
        :return: Type of bus
        :rtype: string
        """
        return self._type_verbose

    ### ----------------------------------------- Setters -----------------------------------------
    @P.setter
    def P(self, p):
        if self._type == 'S': # only settable in S bus
            self._P = p

    @Q.setter
    def Q(self, q):
        if self._type != 'D': # only settable in a S or PV bus
            self._Q = q

    @V.setter
    def V(self, v):
        if self._type == 'D': # only settable in a PQ bus
            self._V = v

    @Th.setter
    def Th(self, th):
        if self._type != 'S': # only settable in a PQ or PV bus
            self._Th = th

    ### ----------------------------------------- Functions -----------------------------------------

    def __init__(self, b: int, p: float, q: float, type: chr, p_gen: int, v_set: float):
        self._b = b - 1 # Bus Number index
        self._type = type # Bus type
        # Create verbose version of type
        if self._type == 'S': self._type_verbose = "Slack"
        elif self._type == 'G': self._type_verbose = "PV"
        else: self._type_verbose = "PQ"

        self._p_load = p / s_base # Load Active Power in pu
        self._q_load = q / s_base # Load Reactive Power in pu
        self._p_gen = p_gen / s_base # Generated Active Power in pu
        self._v_set = v_set # Voltage at the bus in pu

        self._P = None
        if self._type == 'G' or self._type == 'D':
            self._P = self._p_gen - self._p_load

        self._Q = None
        if self._type == 'D':
            self._Q = -self._q_load

        self._V = None
        if self._type == 'S' or self._type == 'G':
            self._V = self._v_set

        self._Th = None
        if self._type == 'S':
            self._Th = 0


    def __repr__(self):
        return (f"{self.__class__.__name__}> Bus: {self._b + 1}, Type: {self.type}")