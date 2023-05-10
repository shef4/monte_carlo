import numpy as np

def metropolis_monte_carlo(ham, conf, T=1, nsweep=1000, nburn=100):
    """
    Perform Metropolis Monte Carlo simulation to obtain thermodynamic properties of a given system.

    Parameters:
    -----------
    ham : object
        An instance of a Hamiltonian class representing the system being studied.
    conf : object
        An instance of a Configuration class representing the initial configuration of the system.
    T : float, optional
        Temperature of the system in units of energy. Default is 1.
    nsweep : int, optional
        The total number of Monte Carlo sweeps to be performed. Default is 1000.
    nburn : int, optional
        The number of sweeps used for thermalization of the system. Default is 100.

    Returns:
    --------
    E_samples : numpy.ndarray
        Array of length nsweep containing the energy samples obtained during the simulation.
    M_samples : numpy.ndarray
        Array of length nsweep containing the magnetization samples obtained during the simulation.
    EE_samples : numpy.ndarray
        Array of length nsweep containing the squared energy samples obtained during the simulation.
    MM_samples : numpy.ndarray
        Array of length nsweep containing the squared magnetization samples obtained during the simulation.
    """
    E_samples = np.zeros(nsweep)
    M_samples = np.zeros(nsweep)
    EE_samples = np.zeros(nsweep)
    MM_samples = np.zeros(nsweep)

    # thermalization
    for _ in range(nburn):
        ham.metropolis_sweep(conf, T=T)

    # accumulation
    E_samples[0] = ham.energy(conf)
    M_samples[0] = conf.get_magnetization()
    MM_samples[0] = M_samples[0] ** 2
    EE_samples[0] = E_samples[0] ** 2

    for si in range(1, nsweep):
        ham.metropolis_sweep(conf, T=T)
        Ei = ham.energy(conf)
        Mi = conf.get_magnetization()

        # update energy and magnetization means
        E_samples[si] = E_samples[si-1] + (Ei - E_samples[si-1]) / (si+1)
        M_samples[si] = M_samples[si-1] + (Mi - M_samples[si-1]) / (si+1)

        # update energy and magnetization squared means
        EE_samples[si] = EE_samples[si-1] + (Ei ** 2 - EE_samples[si-1]) / (si+1)
        MM_samples[si] = MM_samples[si-1] + (Mi ** 2 - MM_samples[si-1]) / (si+1)

    return E_samples, M_samples, EE_samples, MM_samples
