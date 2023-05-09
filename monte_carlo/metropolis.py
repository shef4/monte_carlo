import numpy as np

def metropolis_monte_carlo(ham, conf, T=1, nsweep=1000, nburn=100):
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


