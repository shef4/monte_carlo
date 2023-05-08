Getting Started
===============

Once installed, you can use the package. This example shows how to caclulate the energy and magnetization of a spefic spin configruation.
.. code-block:: python
import monte_carlo

coupling_constant = 1
spin_con = monte_carlo.spin_config("++--+-+++--+-++")
magnet = monte_carlo.magnetization(spin_con)
hamiltonian_energy = monte_carlo.hamiltonian_energy_1d(coupling_constant, spin_con)

