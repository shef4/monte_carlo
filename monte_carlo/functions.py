"""Provide the primary functions."""


def spin_config(spin_string):
    """Converts a string of spin configurations into a list of integers representing the spins.

    Args:
        spin_string (str): A string of spin configurations, where '+' represents spin up and '-' represents spin down.

    Returns:
        list: A list of integers representing the spins, where 1 represents spin up and -1 represents spin down.
    """
    spins = []
    for char in spin_string:
        if char == '+':
            spins.append(1)
        elif char == '-':
            spins.append(-1)
        else:
            print("invalid char, expected '-' or '+'")
            return []
    return spins


def hamiltonian_energy_1d(coupling_const, spins):
    """Calculates the Hamiltonian energy of a 1D spin configuration.

    Args:
        coupling_constant (int): The coupling constant of the Hamiltonian.
        spin_config (list): A list of integers representing the spins of the 1D spin configuration.

    Returns:
        float: The Hamiltonian energy of the 1D spin configuration.
    """
    if len(spins) == 0:
        print ("invalid spin config, list length is 0")
        return 0
    if coupling_const == 0:
        print ("invalid coupling constant, value is 0")
        return 0
    energy = 0
    for i in range(len(spins) - 1):
        energy += spins[i] * spins[i + 1]
    energy *= -coupling_const
    return energy


def magnetization(spins):
    """Calculates the magnetization of a spin configuration.

    Args:
        spin_config (list): A list of integers representing the spins of the spin configuration.

    Returns:
        int: The magnetization of the spin configuration.
    """
    if len(spins) == 0:
        print ("invalid spin config, list length is 0")
        return 0
    return sum(spins)


if __name__ == "__main__":
    print("Please provide spin configuration.\nExample:\n'Spin Configuration: +--+-+--+', '+' is spin up and '-' is spin down.,\n if left empty example value is used as default.")
    spin_str = input("Spin Configuration: ")
    print("Please provide coupling constant.\nExample:\n'Coupling Constant: 1',\n if left empty example value is used as default.")
    coupling_constant = int(input("Coupling Constant: "))
    spin_con = spin_config(spin_str)
    magnet = magnetization(spin_con)
    hamiltonian_energy = hamiltonian_energy_1d(coupling_constant, spin_con)
    print("The energy of the spin configuration",spin_con,"is",hamiltonian_energy,"and it's magnetization is",magnet)
