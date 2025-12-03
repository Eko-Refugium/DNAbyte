import random
from dnabyte.store import SimulateStorage

class Roomtemperature(SimulateStorage):
    # Reference::
    # half-life of DNA: 521 years => lambda = 0.00133
    
    def __init__(self, params, logger=None):
        self.years = params.years

    def simulate(self, data):
        """
        Simulate storage of DNA sequences at room temperature.
        :param assembled_data: An object of class AssembledData.
        :return: A list of stored DNA sequences.
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in data:
            # half-life of DNA: 521 years => lambda = ln(2)/521 ≈ 0.00133 per year
            # decay rate per nucleotide per year
            lambda_decay = 0.00133
            decay_probability = 1 - (1 - lambda_decay)**(len(oligo) * self.years)
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