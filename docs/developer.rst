API Documentation
=================

Requirements
============
   Python 3.x
   pytest (for testing)

Installation
============
Clone the repository or download the code
Install the requirements by running pip install -r requirements-dev.txt

Testing
=======
Run the tests with pytest by running pytest in the project directory.

API documentation
=================

spin_config(spin_string)
========================
Converts a string of spin configurations into a list of integers representing the spins.
Args:
spin_string (str): A string of spin configurations, where '+' represents spin up and '-' represents spin down.
Returns:
list: A list of integers representing the spins, where 1 represents spin up and -1 represents spin down.

hamiltonian_energy_1d(coupling_const, spins)
============================================
Calculates the Hamiltonian energy of a 1D spin configuration.
Args:
coupling_constant (int): The coupling constant of the Hamiltonian.
spin_config (list): A list of integers representing the spins of the 1D spin configuration.
Returns:
float: The Hamiltonian energy of the 1D spin configuration.

magnetization(spins)
====================
Calculates the magnetization of a spin configuration.
Args:
spin_config (list): A list of integers representing the spins of the spin configuration.
Returns:
int: The magnetization of the spin configuration.