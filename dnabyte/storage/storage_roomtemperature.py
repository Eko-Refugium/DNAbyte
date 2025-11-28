import random
from dnabyte.store import SimulateStorage

class Roomtemperature(SimulateStorage):
    # Reference::
    # half-life of DNA: 521 years => lambda = 0.00133
    def simulate(self, assembled_data):
        """
        Simulate storage of DNA sequences at room temperature.
        :param assembled_data: An object of class AssembledData.
        :return: A list of stored DNA sequences.
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in assembled_data:
            decay_probability = 1 - (1 - 5E-3)**(len(oligo) * self.years)
            if random.random() > decay_probability:
                remaining_oligos.append(oligo)
            else:
                strand_breaks += 1

        info = {"number_of_strand_breaks": strand_breaks}
        return remaining_oligos, info
