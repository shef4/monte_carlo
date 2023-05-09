API documentation
=================
Developer Documentation
=======================
BitString class:
The BitString class represents a binary string of a fixed length. It has the following attributes:

N: The length of the binary string.
pbc: A boolean indicating whether to apply periodic boundary conditions when flipping sites.
n_dim: The dimension of the bit string (2^N).
config: A numpy array representing the binary string.
The BitString class provides the following methods:

__init__(self, N=10, pbc=True): Constructs a BitString object with the given length and periodic boundary condition.
__repr__(self): Returns a string representation of the binary string.
__str__(self): Returns a string representation of the binary string.
__getitem__(self, i): Returns the i-th element of the binary string.
initialize(self, M=0): Initializes the binary string with M randomly chosen 1s.
flip_site(self, i): Flips the i-th site of the binary string.
set_int_config(self, int_index): Sets the binary string from an integer index.
get_magnetization(self): Returns the magnetization of the bit string.
set_config(self, conf): Sets the bit string to a given configuration.

metropolis_monte_carlo function:
   metropolis_monte_carlo Function
   The metropolis_monte_carlo function is used for performing Monte Carlo simulations on the Ising model. It has the following parameters:

   ham: An instance of the IsingHamiltonian class.
   conf: An instance of the BitString class representing the initial configuration.
   T: The temperature of the system (default is 1).
   nsweep: The number of Monte Carlo sweeps to perform (default is 1000).
   nburn: The number of thermalization sweeps to perform (default is 100).
   The function returns the following numpy arrays:

   E_samples: The energy samples of the system.
   M_samples: The magnetization samples of the system.
   EE_samples: The energy squared samples of the system.
   MM_samples: The magnetization squared samples of the system.