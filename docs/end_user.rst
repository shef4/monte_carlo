Getting Started
===============

.. _Requirements:
    Python 3.10

.. _Installation:
    Clone the repository or download the code
    Install the requirements by running pip install -r requirements.txt

.. _End-User Documentation:
    CLASSES

.. _BitString class:
    The BitString class represents a binary string of a fixed length. It provides the following methods:

    __init__(self, N=10, pbc=True): Constructs a BitString object with the given length and periodic boundary condition.

    initialize(self, M=0): Initializes the binary string with M randomly chosen 1s.

    flip_site(self, i): Flips the i-th site of the binary string.

    set_int_config(self, int_index): Sets the binary string from an integer index.

    get_magnetization(self): Returns the magnetization of the bit string.

    set_config(self, conf): Sets the bit string to a given configuration. The conf argument should be a list or array of length N.

.. _IsingHamiltonian class:
    The IsingHamiltonian class represents a Hamiltonian of the Ising model. It provides the following methods:

    __init__(self, J=[[()]], mu=np.zeros(1)): Constructs an IsingHamiltonian object with the given coupling matrix and magnetic field.

    energy(self, config): Returns the energy of the Ising model for a given configuration. The config argument should be a BitString object.

.. _metropolis_monte_carlo function:
    The metropolis_monte_carlo function performs Monte Carlo simulations of the Ising model using the Metropolis algorithm. It takes the following arguments:

    ham: An IsingHamiltonian object representing the Hamiltonian of the Ising model.

    conf: A BitString object representing the initial configuration of the system.

    T=1: The temperature of the system (default is 1).

    nsweep=1000: The number of Monte Carlo sweeps to perform (default is 1000).

    nburn=100: The number of sweeps to use for thermalization (default is 100).

    The function returns four arrays of length nsweep: E_samples, M_samples, EE_samples, and MM_samples, representing the energy, magnetization, energy squared, and magnetization squared, respectively, at each Monte Carlo sweep.


.. _Usage:
    To use the code, import the functions and call them with appropriate parameters.

.. code-block:: python

    import monte_carlo
    # Create an instance of the BitString class
    bitstr = monte_carlo.BitString(N=10)

    
    J = [[(1, 1), (2, 1)], [(2, -1)]]
    mu = [0, 0]
    hamiltonian = monte_carlo.IsingHamiltonian(J=J, mu=mu)

    E, M, EE, MM = monte_carlo.metropolis_monte_carlo(hamiltonian, bitstr, T=T, nsweep=8000, nburn=1000)

.. _License:
    This code is released under the MIT License.