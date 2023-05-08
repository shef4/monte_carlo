Getting Started
===============

Requirements
============
Python 3.x
Installation
============
Clone the repository or download the code
Install the requirements by running pip install -r requirements.txt
Usage
=======
To use the code, import the functions and call them with appropriate parameters.

.. code-block:: python

    import monte_carlo
    spin_str = "+--+-+--+"
    coupling_constant = 1

    spin_con = monte_carlo.spin_config(spin_str)
    magnet = monte_carlo.magnetization(spin_con)
    hamiltonian_energy = hamiltonian_energy_1d(coupling_constant, spin_con)
    print("The energy of the spin configuration",spin_con,"is",hamiltonian_energy,"and it's magnetization is",magnet)


License
=======
This code is released under the MIT License.