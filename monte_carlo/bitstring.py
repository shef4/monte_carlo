import numpy as np
import random

class BitString:
    """
    A class representing a binary string of length N.

    Attributes:
    - N (int): The length of the bitstring.
    - pbc (bool): Whether the bitstring should be treated as a periodic boundary condition.
    - n_dim (int): The total number of possible bitstrings of length N.
    - config (numpy.ndarray): The binary array representing the current state of the bitstring.

    Methods:
    - __init__(self, N=10, pbc=True): Initializes a new BitString object.
    - __repr__(self): Returns a string representation of the BitString object.
    - __str__(self): Returns a string representation of the binary configuration of the bitstring.
    - __getitem__(self, i): Returns the binary value at index i in the bitstring.
    - initialize(self, M=0): Randomly initializes the bitstring with M non-zero values.
    - flip_site(self, i): Flips the binary value at index i in the bitstring.
    - set_int_config(self, int_index): Sets the configuration of the bitstring to the binary representation of int_index.
    - get_magnetization(self): Calculates and returns the magnetization of the bitstring.
    - set_config(self, conf): Sets the configuration of the bitstring to conf.
    """
    def __init__(self, N=10, pbc=True):
        """
        Constructs a BitString object with the given length and periodic boundary condition.

        Parameters
        ----------
        N : int, optional
            The length of the binary string (default is 10).
        pbc : bool, optional
            Whether to apply periodic boundary conditions when flipping sites (default is True).
        """
        self.N = N
        self.pbc = pbc
        self.n_dim = 2 ** self.N
        self.config = np.zeros(N, dtype=int)

    def __repr__(self):
        """
        Returns a string representation of the binary string.
        """
        return f"BitString({self.config})"

    def __str__(self):
        """
        Returns a string representation of the binary string.
        """
        return "".join(str(e) for e in self.config)

    def __getitem__(self, i):
        """
        Returns the i-th element of the binary string.

        Parameters
        ----------
        i : int
            The index of the element to get.

        Returns
        -------
        int
            The value of the i-th element (0 or 1).
        """
        return self.config[i]

    def initialize(self, M=0):
        """
        Initializes the binary string with M randomly chosen 1s.

        Parameters
        ----------
        M : int, optional
            The number of 1s to initialize the binary string with (default is 0).
        """
        self.config = np.zeros(self.N, dtype=int)
        random_indices = random.sample(range(0, self.N), M)
        for i in random_indices:
            self.config[i] = 1

    def flip_site(self, i):
        """
        Flips the i-th site of the binary string.

        Parameters
        ----------
        i : int
            The index of the site to flip.
        """
        self.config[i] = (self.config[i]+1)%2

    def set_int_config(self, int_index):
        """
        Sets the binary string from an integer index.

        Parameters
        ----------
        int_index : int
            The integer index to set the binary string from.
        """
        self.config = np.array([int(i) for i in np.binary_repr(int_index, width=self.N)])

    def get_magnetization(self):
        """
        Returns the magnetization of the bit string.
        """
        return np.sum(2 * self.config - 1)

    def set_config(self, conf):
        """
        Sets the bit string to a given configuration.
        
        Parameters:
        -----------
        conf: list or array
            The configuration to set the bit string to. Must be of length N.
        """
        assert len(conf) == self.N
        self.config = np.array(conf)
