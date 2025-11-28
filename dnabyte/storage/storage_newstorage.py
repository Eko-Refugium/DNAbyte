from dnabyte.store import SimulateStorage
from dnabyte.data_classes.insilicodna import InSilicoDNA
import random

class Newstorage(SimulateStorage):
    def simulate(self, assembled_data):
        """
        Simulate storage of DNA sequences in cryogenic conditions.
        :param assembled_data: An object of class AssembledData.
        :return: A tuple (stored_data, strand_breaks).
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in assembled_data:
            decay_probability = 1 - (1 - 1E-8)**(len(oligo) * self.years)
            if random.random() > decay_probability:
                remaining_oligos.append(oligo)
            else:
                strand_breaks += 1

        stored_data = InSilicoDNA(assembled_data=remaining_oligos)
        info = {"number_of_strand_breaks": strand_breaks}
        return stored_data, info