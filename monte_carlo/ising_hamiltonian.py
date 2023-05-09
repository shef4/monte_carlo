import numpy as np
import random

class IsingHamiltonian:
    def __init__(self, J=[[()]], mu=np.zeros(1)):
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
        delta_e = 0.0
        delta_si = 2
        if config[i] == 1:
            delta_si = -2

        for j, j_val in enumerate(self.nodes[i]):
            delta_e += (2.0*config[j_val]-1.0) * self.js[i][j] * delta_si

        delta_e += self.mu[i] * delta_si
        return delta_e

    def metropolis_sweep(self, conf, T=1.0):
        for site_i in range(conf.N):
            delta_e = self.delta_e_for_flip(site_i, conf.config)

            accept = True
            if delta_e > 0.0:
                rand_comp = random.random()
                if rand_comp > np.exp(-delta_e/T):
                    accept = False

            if accept:
                conf.config[site_i] = 1 - conf.config[site_i]

        return conf

    def compute_average_values(self, conf, T):
        energy, magnetization, partition, energy_squared, magnetization_squared = 0.0, 0.0, 0.0, 0.0, 0.0

        for i in range(conf.n_dim):
            conf.set_int_config(i)
            energy_i = self.energy(conf)
            partition_i = np.exp(-energy_i/T)
            energy += energy_i * partition_i
            energy_squared += energy_i * energy_i * partition_i
            magnetization_i = np.sum(2*conf.config-1)
            magnetization += magnetization_i * partition_i
            magnetization_squared += magnetization_i * magnetization_i * partition_i
            partition += partition_i

        energy /= partition
        magnetization /= partition
        energy_squared /= partition
        magnetization_squared /= partition

        heat_capacity = (energy_squared - energy*energy)/(T*T)
        magnetic_susceptibility = (magnetization_squared - magnetization*magnetization)/T
        return energy, magnetization, heat_capacity, magnetic_susceptibility

