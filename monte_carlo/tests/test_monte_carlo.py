"""
Unit and regression test for the monte_carlo package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import monte_carlo


def test_monte_carlo_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "monte_carlo" in sys.modules
    spin_list = monte_carlo.spin_config("+-++---+")
    assert (spin_list == [1,-1,1,1,-1,-1,-1,1])
    spin_list = monte_carlo.spin_config("+-+Gkc--+")
    assert (not spin_list)

    energy = monte_carlo.hamiltonian_energy_1d(1, [1,-1,1,1,-1,-1,-1,1])
    assert (energy == 1)
    energy = monte_carlo.hamiltonian_energy_1d(0, [1,-1,1,1,-1,-1,-1,1])
    assert (energy == 0)
    energy = monte_carlo.hamiltonian_energy_1d(0, [])
    assert (energy == 0)
    
    magnet = monte_carlo.magnetization([1,-1,1,1,-1,1,-1,1])
    assert (magnet == 2)

