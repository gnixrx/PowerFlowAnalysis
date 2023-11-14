import time
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
        :return: Time and number of iterations it took to complete the analysis.
        :rtype: string
        """
        start = time.perf_counter()
        self.solve_implicit()
        self.solve_explicit()
        end = time.perf_counter()

        return f"Update took {end - start:0.5f}s in {self._iterations} iterations."

    def solve_implicit(self):
        """
        Solves the implicit equations for power analysis.
        :return: None
        :rtype: None
        """
        # Start implicit equations
        convergence = 1
        self._iterations = 0
        while convergence >= 0.001:  # When the minimum is over the floor accuracy
            mm_vector = self.calc_mismatch()

            # Build inverse Jacobian matrix
            J = np.bmat([[self.H(), self.M()], [self.N(), self.L()]])
            J_inv = -1 * np.linalg.inv(J)

            # Find delta theta, delta V
            delta_Th_V = np.ravel(J_inv @ mm_vector)

            # Update theta, and V by the delta
            i = 0
            for bus in self._pq_pv:
                bus.Th = bus.Th + delta_Th_V[i]
                if bus._type == 'D':
                    bus.V = bus.V + delta_Th_V[i + self._pq_pv.size]
                i += 1

            convergence = np.max(np.abs(delta_Th_V))  # Check to see if need for more adjustments
            self._iterations += 1

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
        V_k = k.V
        V_i = i.V
        Th_ki = k.Th - i.Th
        G_ki = self._G_mat[k._b, i._b]
        B_ki = self._B_mat[k._b, i._b]
        return V_k, V_i, Th_ki, G_ki, B_ki


    def H(self):
        """
        Creates the H matrix the upper left of the Jacobian matrix.
        :return: H matrix
        :rtype: Array
        """
        H = np.zeros((self._pq_pv.size, self._pq_pv.size))

        k = 0
        for bus_k in self._pq_pv:
            H_kk = []
            # Iterate through connections to k diagonal
            for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                bus_g = self._bus_data[self._Y_col[g][0]]
                V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_g)

                if bus_k._b != bus_g._b:
                    H_kk.append(V_k * V_i * (-G_ki * np.sin(Th_ki) + B_ki * np.cos(Th_ki)))
            H[k, k] = np.sum(H_kk)

            # Iterate connections between ki nondiagonal
            i = 0
            for bus_i in self._pq_pv:
                # bus_i = self._pq_pv[i]
                V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    H[k, i] = V_k * V_i * (G_ki * np.sin(Th_ki) - B_ki * np.cos(Th_ki))

                i += 1
            k += 1

        return H

    def M(self):
        """
        Creates the M matrix the upper right of the Jacobian matrix.
        :return: M matrix
        :rtype: Array
        """
        M = np.zeros((self._pq_pv.size, self._pq.size))

        k = 0
        for bus_k in self._pq_pv:
            # Iterate through connections to kk
            if k < self._pq.size:
                M_kk = []
                M_kk_front = 0
                for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                    bus_g = self._bus_data[self._Y_col[g][0]]
                    V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_g)

                    if bus_k._b != bus_g._b:
                        M_kk.append(V_i * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki)))
                    else:
                        M_kk_front = 2 * G_ki * V_k
                M[k, k] = M_kk_front + np.sum(M_kk)

            # Iterate connections between ki
            i = 0
            for bus_i in self._pq:
                V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    M[k, i] = V_k * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki))
                i += 1
            k += 1

        return M

    def N(self):
        """
        Creates the N matrix the lower left of the Jacobian matrix.
        :return: N matrix
        :rtype: Array
        """
        N = np.zeros((self._pq.size, self._pq_pv.size))

        k = 0
        for bus_k in self._pq:
            N_kk = []

            # Iterate through connections to kk
            for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                bus_g = self._bus_data[self._Y_col[g][0]]
                V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_g)

                if bus_k._b != bus_g._b:
                    N_kk.append(V_k * V_i * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki)))
            N[k, k] = np.sum(N_kk)

            i = 0
            # Iterate connections between ki
            for bus_i in self._pq_pv:
                V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    N[k, i] = V_k * V_i * (-G_ki * np.cos(Th_ki) - B_ki * np.sin(Th_ki))
                i += 1
            k += 1

        return N

    def L(self):
        """
        Creates the L matrix the lower right of the Jacobian matrix.
        :return: L matrix
        :rtype: Array
        """
        L = np.zeros((self._pq.size, self._pq.size))

        k = 0
        for bus_k in self._pq:
            L_kk = []
            L_kk_front = 0
            # Iterate through connections to kk
            for g in np.argwhere(self._Y_row == bus_k._b):  # Each connection where admittance != 0)
                bus_g = self._bus_data[self._Y_col[g][0]]
                V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_g)

                if bus_k._b != bus_g._b:
                    L_kk.append(V_i * (G_ki * np.sin(Th_ki) - B_ki * np.cos(Th_ki)))
                else:
                    L_kk_front = -2 * B_ki * V_k
            L[k, k] = L_kk_front + np.sum(L_kk)

            # Iterate connections between ki
            i = 0
            for bus_i in self._pq:
                V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus_k, bus_i)

                if bus_k._b != bus_i._b:
                    L[k, i] = V_k * (G_ki * np.sin(Th_ki) - B_ki * np.cos(Th_ki))

                i += 1
            k += 1

        return L


    def solve_explicit(self):
        """
        Solves the explicit equations with knowns from implicit equations.
        :return: None
        :rtype: None
        """
        for bus in self._s_pv:
            if bus._type == 'S':
               bus.P = np.sum(self.power_flow_eq(True, bus))
            bus.Q = np.sum(self.power_flow_eq(False, bus))

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
            V_k, V_i, Th_ki, G_ki, B_ki = self.values_ki(bus, bus_i)

            if isP:
                Pk_Qk.append(V_k * V_i * (G_ki * np.cos(Th_ki) + B_ki * np.sin(Th_ki)))
            else: # isQ
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
            if bus._type == 'D':
                PQ[P_i + self._pq_pv.size] = np.sum(self.power_flow_eq(False, bus)) - bus._Q
            P_i += 1 # Next PQ/PV bus

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
        # Ordered bus_data
        self._s_pv = np.array([bus for bus in bus_data if bus._type == 'S' or bus._type == 'G'])
        self._pq = np.array([bus for bus in bus_data if bus._type == 'D'])
        self._pq_pv = np.concatenate((self._pq, np.array([bus for bus in bus_data if bus._type == 'G'])))

    def __repr__(self):
        return f"{self.__class__.__name__}> Number of Nodes: {self._node_max}"