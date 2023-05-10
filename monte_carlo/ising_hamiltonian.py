import numpy as np
import random

class IsingHamiltonian:
    """
    A class representing an Ising Hamiltonian system with the following attributes and methods:

    Attributes:

    J: A 2D array of size (N x N) representing the interaction strength between spins, where N is the number of spins. J[i][j] represents the interaction strength between spins i and j.
    mu: A 1D array of size N representing the external field strength for each spin.
    nodes: A list of length N, where each element is a 1D array containing the indices of the spins that interact with the corresponding spin.
    js: A list of length N, where each element is a 1D array containing the corresponding interaction strengths for the spins in the corresponding node.
    Methods:

    init(self, J=[[()]], mu=np.zeros(1)): Constructs an instance of the IsingHamiltonian class with given interaction strength J and external field strength mu. Initializes nodes and js attributes.
    energy(self, config): Computes the energy of the system for a given configuration of spins.
    delta_e_for_flip(self, i, config): Computes the change in energy due to a flip of the spin at index i for a given configuration of spins.
    metropolis_sweep(self, conf, T=1.0): Performs a single Metropolis sweep of the Ising Hamiltonian system for a given configuration of spins and temperature T.
    compute_average_values(self, conf, T): Computes the average energy, magnetization, specific heat, and magnetic susceptibility of the Ising Hamiltonian system for a given configuration of spins and temperature T.
    """
    def __init__(self, J=[[(0)]], mu=np.zeros(1)):
        """
        Initializes the IsingHamiltonian object with the given coupling coefficients and external magnetic field values.

        Parameters:
        J (list of lists of tuples): The coupling coefficients between each pair of spins in the Ising model.
        mu (ndarray): The values of the external magnetic field at each spin site.

        Returns:
        None
        """
        self.J = J
        self.mu = mu

        self.nodes = []
        self.js = []

        for i in range(len(self.J)):
            node_arr = np.zeros(len(self.J[i]), dtype=int)
            j_arr = np.zeros(len(self.J[i]))
            for j, (node, j_val) in enumerate(self.J[i]):
                node_arr[j] = node
                j_arr[j] = j_val
            self.nodes.append(node_arr)
            self.js.append(j_arr)

        self.mu = np.array(self.mu)

    def energy(self, config):
        """
        Calculates the energy of the given spin configuration based on the current Ising Hamiltonian.

        Parameters:
        config (IsingConfig): The spin configuration for which the energy needs to be calculated.

        Returns:
        float: The energy of the spin configuration.
        """
        if len(config.config) != len(self.J):
            raise ValueError("The configuration's length does not match the Hamiltonian's dimension.")

        energy = 0.0
        for i in range(config.N):
            for j, (node, j_val) in enumerate(self.J[i]):
                if node < i:
                    continue
                if config[i] == config[node]:
                    energy += j_val
                else:
                    energy -= j_val

        energy += np.dot(self.mu, 2*config.config-1)

        return energy

    def delta_e_for_flip(self, i, config):
        """
        Calculates the change in energy if the i-th spin in the configuration is flipped.

        Parameters:
        i (int): The index of the spin in the configuration that needs to be flipped.
        config (IsingConfig): The spin configuration for which the energy change needs to be calculated.

        Returns:
        float: The change in energy due to flipping the i-th spin.
        """
        delta_e = 0.0
        delta_si = 2
        if config[i] == 1:
            delta_si = -2

        for j, j_val in enumerate(self.nodes[i]):
            delta_e += (2.0*config[j_val]-1.0) * self.js[i][j] * delta_si

        delta_e += self.mu[i] * delta_si
        return delta_e

    def metropolis_sweep(self, conf, T=1.0):
        """
        Executes a single Metropolis sweep of the Ising model at the given temperature.

        Parameters:
        conf (IsingConfig): The initial spin configuration for the Metropolis sweep.
        T (float): The temperature at which the Metropolis sweep is to be performed.

        Returns:
        IsingConfig: The spin configuration after the Metropolis sweep.
        """
        for site_i in range(conf.N):
            delta_e = self.delta_e_for_flip(site_i, conf.config)      
            accept = True
            if delta_e > 0.0:
                # prob_trans = np.exp(-delta_e/T)
                rand_comp = random.random()
                if rand_comp > np.exp(-delta_e/T):
                    accept = False
            if accept:
                if conf.config[site_i] == 0:
                    conf.config[site_i] = 1
                else:
                    conf.config[site_i] = 0
        return conf

    
    def compute_average_values(self, conf, T):
        """
        Computes the average energy, magnetization, specific heat, and magnetic susceptibility of the Ising model at the given temperature.

        Parameters:
        conf (IsingConfig): The spin configuration for which the average values need to be calculated.
        T (float): The temperature at which the average values need to be calculated.

        Returns:
        tuple: A tuple containing the average energy, magnetization, specific heat, and magnetic susceptibility of the Ising model.
        """
        E  = 0.0
        M  = 0.0
        Z  = 0.0
        EE = 0.0
        MM = 0.0

        for i in range(conf.n_dim):
            conf.set_int_config(i)
            Ei = self.energy(conf)
            Zi = np.exp(-Ei/T)
            E += Ei*Zi
            EE += Ei*Ei*Zi
            Mi = np.sum(2*conf.config-1)
            M += Mi*Zi
            MM += Mi*Mi*Zi
            Z += Zi
        
        E = E/Z
        M = M/Z
        EE = EE/Z
        MM = MM/Z
        
        HC = (EE - E*E)/(T*T)
        MS = (MM - M*M)/T
        return E, M, HC, MS

