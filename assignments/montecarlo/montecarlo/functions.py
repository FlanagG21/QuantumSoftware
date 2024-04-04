"""Provide the primary functions."""
k =  1
def canvas(with_attribution=True):
    """
    Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """

    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())

class BitString:
    """
    Simple class to implement a config of bits
    """
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int)

    def __repr__(self):
        pass

    def __eq__(self, other):
        if other is None :
            return False
        if type(self) != type(other) :
            return False
        if self.N != other.N :
            return False              
        return np.array_equal(self.config, other.config)


    
    def __len__(self):
        return self.N

    def on(self):
        return self.config.sum()

    
    def off(self):
        return self.N - self.config.sum()
    
    def flip_site(self,i):
        self.config[i] = (self.config[i] + 1) % 2 
    
    def int(self):
        num = 0
        for i, bit in enumerate(reversed(self.config)):
            num += bit * (2 ** i)
        return num

    def set_config(self, s:list[int]):
        for i, bit in enumerate(s):
            self.config[i + (self.N - len(s))] = bit
        
    def set_int_config(self, dec:int):        
        for i in range(len(self.config)):
            self.config[i] = dec % 2
            dec //= 2
        self.config = self.config[::-1]
        



    def __str__(self):
        builder = ""
        for bit in self.config:
            builder += str(bit)
        return builder


def energy(bs: BitString, G: nx.Graph):
    """Compute energy of configuration, `bs`

        .. math::
            E = \\left<\\hat{H}\\right>

    Parameters
    ----------
    bs   : Bitstring
        input configuration
    G    : Graph
        input graph defining the Hamiltonian
    Returns
    -------
    energy  : float
        Energy of the input configuration
"""
def energy(bs: BitString, G: nx.Graph):
    energy = 0.0
    for u,v in G.edges():
        if bs.config[u] == bs.config[v]:
            energy += G.edges[u,v]["weight"]
        else:
            energy -= G.edges[u,v]["weight"]
    return energy



def compute_average_values(bs:BitString, G: nx.Graph, T: float):
    """
    Compute the average value of Energy, Magnetization, 
    Heat Capacity, and Magnetic Susceptibility 

        .. math::
            E = \\left<\\hat{H}\\right>

    Parameters
    ----------
    bs   : Bitstring
        input configuration
    G    : Graph
        input graph defining the Hamiltonian
    T    : float
        temperature of the system
    Returns
    -------
    energy  : float
    magnetization  : float
    heat capacity  : float
    magnetic susceptibility  : float
    """
    E = 0
    M = 0
    HC = 0
    MS = 0

    z, z2, zm, zm2 = boltzmanDenominator(bs, G, T)
    for i in range(2 ** bs.N):
        bs.set_int_config(i)
        currEnergy = energy(bs, G)
        p = np.e ** ((-1 / (k * T))* currEnergy)
        p /= z
        E += currEnergy*p
        currEnergySquared = currEnergy ** 2
        HC += currEnergySquared * p
        mag = bs.on() - bs.off()
        M += mag * p
        magSQR = mag**2
        MS += magSQR * p
    MS = (MS - M**2) * T ** -1        
    HC = (HC - E**2) * T ** -2
    
    return E, M, HC, MS

def boltzmanDenominator(bs:BitString, G: nx.Graph, T: float):
    z = 0
    z2 = 0
    zm = 0
    zm2 = 0
    for i in range(2 ** bs.N):
        bs.set_int_config(i)
        currEnergy = energy(bs, G)
        z += np.e ** ((-1 / (k * T))* currEnergy)
        currEnergySquared = currEnergy ** 2
        z2 += np.e ** ((-1 / (k * T)) * currEnergySquared)
        mag = bs.on() - bs.off()
        magSQR = mag ** 2
        zm += np.e ** ((-1 / (k * T)) * mag)
        zm2 += np.e ** ((-1 / (k * T)) * magSQR)
    currEnergy = energy(bs, G)
    z += np.e ** ((-1 / (k * T)) * currEnergy)
    print((-1 / (k * T)))
    return z, z2, zm, zm2