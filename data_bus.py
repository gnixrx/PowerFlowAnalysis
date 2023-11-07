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

TODO:

"""

class BusData:

    def __init__(self, b: int, p: float, q: float, type: chr, p_gen: int, v_set: float):
        self._b = b - 1 # Bus Number index
        self._type = type # Bus type
        # Create verbose version of type
        if self._type == 'S': self._type_verbose = "Slack bus"
        elif self._type == 'G': self._type_verbose = "PV bus"
        else: self._type_verbose = "PQ bus"

        self._p_load = p / s_base # Load Active Power in pu
        self._q_load = q / s_base # Load Reactive Power in pu
        self._p_gen = p_gen / s_base # Generated Active Power in pu

        if self._type == 'S' or self._type == 'PV':
            self._v_set = v_set # Reference voltage at the bus in pu
        else:
            self._v_set = None


    def __repr__(self):
        return f"{self.__class__.__name__}> Bus Num: {self._b + 1}, Type: {self._type_verbose}, "