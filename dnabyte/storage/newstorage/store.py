from dnabyte.store import SimulateStorage
import random

class Newstorage(SimulateStorage):

    def __init__(self, params, logger=None):
        self.years = params.years

    def simulate(self, data):
        """
        Simulate storage of DNA sequences in cryogenic conditions.
        :param data: An object of class InSilicoDNA.
        :return: A tuple (stored_data, strand_breaks).
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in data:
            decay_probability = 1 - (1 - 1E-8)**(len(oligo) * self.years)
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