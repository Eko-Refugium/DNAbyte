import random
from dnabyte import InSilicoDNA, params
from dnabyte.store import SimulateStorage

class Biogene(SimulateStorage):
    # biogene: (anhydrous and anoxic atmosphere maintained inside hermetic capsules)
    # 1E-7/nt/yr =?= 1 cut/century/100 000 nucleotides = 
    # Reference: 
    # Coudy, Delphine, et al. "Long term conservation of DNA at ambient temperature. 
    # Implications for DNA data storage." PLoS One 16.11 (2021): e0259868.

    def __init__(self, params, logger=None):
        self.years = params.years

    def simulate(self, data):
        """
        Simulate storage of DNA sequences in a biogene storage medium.
        :param assembled_data: An object of class AssembledData.
        :return: A list of stored DNA sequences.
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in data:
            decay_probability = 1 - (1 - 1E-7)**(len(oligo) * self.years)
            if random.random() > decay_probability:
                remaining_oligos.append(oligo)
            else:
                strand_breaks += 1

        info = {"number_of_strand_breaks": strand_breaks}

        return remaining_oligos, info
    
def attributes(params):
    if 'years' not in params.__dict__ or params.years is None:
        years = 100
    else:
        years = params.years

    return {"years": years}