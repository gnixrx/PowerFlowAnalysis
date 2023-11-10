import numpy as np
from scipy.sparse import csr_array
from data_bus import BusData
"""
The power analysis class interperates the inputs from the node/bus and line data and uses the power analysis equations
to find the power injected at the node/bus.
First the class solves the mismatch equations for PQ and PV busses:
    Pk = sum from i=1 to N ( Vk Vi (Gki cos Thki + Bki sin Thki)
    Qk = sum from i=1 to N ( Vk Vi (Gki cos Thki - Bki sin Thki)
:param bus_data: A list of BusData which represent the  nodes in the graph of the network in the power system.
:type bus_data: np.Array of BusData
:param y_matrix: Sparse admittance matrix
:type y_matrix: csr_array
"""

class PowerAnalysis:
    def update(self):
        """
        Run the power system analysis.
        :return: None
        :rtype: None
        """
        self.solve_implicit()
        self.solve_explicit()

    def solve_implicit(self):
        """
        Solves the implicit equations for power analysis.
        :return: None
        :rtype: None
        """
        # Start implicit equations
        mm_vector = self.calc_mismatch().reshape(-1, 1)
        accuracy = np.max(np.abs(mm_vector))  # Initialize a value to start
        while accuracy >= 0.001:  # When the minimum is over the floor accuracy

            # Build inverse Jacobian matrix
            J = np.bmat([[self.H(), self.M()], [self.N(), self.L()]])
            J_inv = -1 * np.linalg.inv(J)

            # Find delta theta, delta V
            delta_Th_V = np.ravel(J_inv @ mm_vector)
            # Update theta, and V by the delta
            for i in range(self._pq_pv.size):
                bus_i = self._pq_pv[i]
                if bus_i._type == 'D':
                    bus_i._V = bus_i._V + delta_Th_V[i + self._pq.size]
                bus_i._Th = bus_i._Th + delta_Th_V[i]

            mm_vector = self.calc_mismatch().reshape(-1, 1)
            accuracy = np.max(np.abs(mm_vector))  # Check to see if need for more adjustments


    def values_ki(self, k, i):
        """
        Return values for Voltage, Voltage Angle, and Acceptance
        :param k: Bus k
        :type k: BusData
        :param i: Bus i
        :type i: Bus Data
        :return: Voltage_k, Voltage_i, Conductance_ki, Susceptance_ki, Theta_delta_ki
        :rtype: float, float, float, float, float
        """
        V_k = k._V
        V_i = i._V
        B_ki = self._B_mat[k._b, i._b]
        G_ki = self._G_mat[k._b, i._b]
        Th_ki = k._Th - i._Th
        return V_k, V_i, B_ki, G_ki, Th_ki


    def H(self):
        H = np.zeros((self._pq_pv.size, self._pq_pv.size))

        for k in range(self._pq_pv.size):
            H_kk = []
            bus_k = self._pq_pv[k]

            # Iterate connections between ki
            for i in range(self._pq_pv.size):
                bus_i = self._pq_pv[i]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    H[k, i] = V_k * V_i * (G_ki * np.sin(Th_ki) - B_ki * np.cos(Th_ki))

            # Iterate through connections to kk
            for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                bus_g = self._bus_data[self._Y_col[g][0]]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_g)

                if bus_k._b != bus_g._b:
                    H_kk.append(V_k * V_i * (-G_ki * np.sin(Th_ki) + B_ki * np.cos(Th_ki)))
                H[k, k] = np.sum(H_kk)

        return H

    def M(self):
        M = np.zeros((self._pq_pv.size, self._pq.size))

        for k in range(self._pq_pv.size):
            bus_k = self._pq_pv[k]
            M_kk = []
            M_kk_front = 0
            # Iterate connections between ki
            for i in range(self._pq.size):
                bus_i = self._pq_pv[i]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    M[k, i] = V_k * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki))

            # Escape early because pq_pv.size > pq.size
            if k >= self._pq.size:
                return M

            # Iterate through connections to kk
            for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                bus_g = self._bus_data[self._Y_col[g][0]]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_g)

                if bus_k._b != bus_g._b:
                    M_kk.append(V_i * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki)))
                else:
                    M_kk_front = 2 * G_ki * V_k
                M[k, k] = M_kk_front + np.sum(M_kk)

        return M

    def N(self):
        N = np.zeros((self._pq.size, self._pq_pv.size))

        for k in range(self._pq.size):
            N_kk = []
            bus_k = self._pq_pv[k]

            # Iterate connections between ki
            for i in range(self._pq_pv.size):
                bus_i = self._pq_pv[i]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    N[k, i] = V_k * V_i * (-G_ki * np.cos(Th_ki) - B_ki * np.sin(Th_ki))

            # Iterate through connections to kk
            for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                bus_g = self._bus_data[self._Y_col[g][0]]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_g)

                if bus_k._b != bus_g._b:
                    N_kk.append(V_k * V_i * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki)))
                N[k, k] = np.sum(N_kk)

        return N

    def L(self):
        L = np.zeros((self._pq.size, self._pq.size))

        for k in range(self._pq.size):
            bus_k = self._pq_pv[k]
            L_kk = []
            L_kk_front = 0
            # Iterate connections between ki
            for i in range(self._pq.size):
                bus_i = self._pq_pv[i]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    L[k, i] = V_k * (G_ki * np.sin(Th_ki) - B_ki * np.cos(Th_ki))

            # Iterate through connections to kk
            for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                bus_g = self._bus_data[self._Y_col[g][0]]
                V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus_k, bus_g)

                if bus_k._b != bus_g._b:
                    L_kk.append(V_i * (G_ki * np.sin(Th_ki) - B_ki * np.cos(Th_ki)))
                else:
                    L_kk_front = -2 * B_ki * V_k
                L[k, k] = L_kk_front + np.sum(L_kk)

        return L


    def solve_explicit(self):
        """
        Solves the explicit equations with knowns from implicit equations.
        :return: None
        :rtype: None
        """
        for i in range(self._s_pv.size):
            bus_i = self._s_pv[i]

            if bus_i._type == 'S':
               bus_i._P = np.sum(self.power_flow_eq(True, bus_i))
            bus_i._Q = np.sum(self.power_flow_eq(False, bus_i))

    def power_flow_eq(self, isP, bus):
        """
        Calculates the mismatch equations for a bus.
        :param isP: Calcluate the P equations if this is true, else the Q equations
        :type isP: boolean
        :param bus: The bus which you want the equations for.
        :type bus: BusData
        :return: List of mismatch equation values to sum for the bus
        :rtype: List(float)
        """
        Pk_Qk = []
        for col_i in np.argwhere(self._Y_row == bus._b):  # Each connection where admittance != 0)
            bus_i = self._bus_data[self._Y_col[col_i][0]]
            V_k, V_i, B_ki, G_ki, Th_ki = self.values_ki(bus, bus_i)

            if isP:
                Pk_Qk.append(V_k * V_i * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki)))
            else:
                Pk_Qk.append(V_k * V_i * (G_ki * np.sin(Th_ki) - B_ki * np.cos(Th_ki)))

        return Pk_Qk

    def calc_mismatch(self):
        """
        Calculates the values for mismatch equations in the power system given the conditions of the Bus Data.
        :return: The values for mismatch equations in a vector
        :rtype: np.array(float)
        """
        # for each P and Q equation or PQ and PV bus calculate the resulting P and Q mismatch numbers
        PQ = np.empty(self._pq_pv.size + self._pq.size)

        # Iterate over the PQ bus
        P_i = 0
        for bus in self._pq_pv:
            # Calculate Active Power for PQ/PV bus
            PQ[P_i] = np.sum(self.power_flow_eq(True, bus)) - bus._P

            P_i += 1 # Next PQ/PV bus

        for bus in self._pq: # Each pq bus
            # Calculate Reactive Power for PQ Bus
            PQ[P_i] = np.sum(self.power_flow_eq(False, bus)) - bus._Q

            P_i += 1

        return PQ


    def get_node_max(self):
        """
        Returns the maximum node numbers in the network from the bus data
        :return: Maximum node
        :rtype: int
        """
        return self._node_max

    def __init__(self, bus_data: np.array, y_matrix: csr_array):
        # Create node count
        self._node_max = max([bus._b for bus in bus_data]) + 1

        # Set up Y matrix
        self._Y_row, self._Y_col = y_matrix.nonzero()
        self._G_mat, self._B_mat = np.real(y_matrix), np.imag(y_matrix)

        # Set up bus_data
        self._bus_data = bus_data
        self._s_pv = np.array([bus for bus in bus_data if bus._type == 'S' or bus._type == 'G'])
        self._pq = np.array([bus for bus in bus_data if bus._type == 'D'])
        self._pq_pv = np.array([bus for bus in bus_data if bus._type == 'D' or bus._type == 'G'])

    def __repr__(self):
        return f"{self.__class__.__name__}> Number of Nodes: {self._node_max}"