"""
Unit and regression test for the monte_carlo package.
"""

# Import package, test suite, and other packages as needed
import numpy as np
import random
import monte_carlo



def test_bitstring_class():
    # Create an instance of the BitString class
    bitstr = monte_carlo.BitString(N=10)

    # Test the __repr__ method
    assert repr(bitstr) == "BitString([0 0 0 0 0 0 0 0 0 0])"

    # Test the __str__ method
    assert str(bitstr) == "0000000000"

    # Test the __getitem__ method
    assert bitstr[0] == 0

    # Test the initialize method
    bitstr.initialize(M=5)
    assert np.sum(bitstr.config) == 5

    # Test the flip_site method
    bitstr.flip_site(0)
    assert bitstr[0] == 1

    # Test the set_int_config method
    bitstr.set_int_config(10)
    value = np.binary_repr(10, width= bitstr.N)
    for i in range(bitstr.N):
        assert bitstr.config[i] == int(value[i])

    # Test the get_magnetization method
    assert bitstr.get_magnetization() == -6

    # Test the set_config method
    value = [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
    bitstr.set_config(value)
    for i in range(bitstr.N):
        assert bitstr.config[i] == int(value[i])


def test_IsingHamiltonian():
    # Initialize Hamiltonian
    J = [[(1, 1), (2, 1)], [(2, -1)]]
    mu = [0, 0]
    hamiltonian = monte_carlo.IsingHamiltonian(J=J, mu=mu)

    # Test delta_e_for_flip method
    config = monte_carlo.BitString(N=3)
    config.initialize(M=2)
    delta_e = hamiltonian.delta_e_for_flip(0, config)
    expected_delta_e = 0
    assert delta_e == expected_delta_e


def test_metropolis():
    # Set random seed for reproducibility
    random.seed(2)

    # Define system parameters
    N = 20
    T = 2
    Jval = 1.0
    mu = [.1 for i in range(N)]
    J = []
    for i in range(N):
        J.append([((i+1) % N, Jval), ((i-1) % N, Jval)])
    ham = monte_carlo.IsingHamiltonian(J=J, mu=mu)

    # Test energy conservation
    conf = monte_carlo.BitString(N=N)
    E, _, _, _ = monte_carlo.metropolis_monte_carlo(ham, conf, T=T, nsweep=8000, nburn=1000)

    # Test that the final energy is within expected range
    assert np.isclose(-9.31, E[-1], rtol=1e-2)

    # Test magnetization conservation
    conf = monte_carlo.BitString(N=N)
    _, M, _, _ = monte_carlo.metropolis_monte_carlo(ham, conf, T=T, nsweep=8000, nburn=1000)

    # Test specific heat calculation
    _, _, EE, _ = monte_carlo.metropolis_monte_carlo(ham, conf, T=T, nsweep=8000, nburn=1000)
    HC = (EE[-1] - E[-1]**2)/T**2

    # Test magnetic susceptibility calculation
    _, _, _, MM = monte_carlo.metropolis_monte_carlo(ham, conf, T=T, nsweep=8000, nburn=1000)
    MS = (MM[-1] - M[-1]**2)/T