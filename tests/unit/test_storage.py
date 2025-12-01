import unittest
from dnabyte.oligo import Oligo
from backup.data import AssembledData, StoredData
from dnabyte.store import SimulateStorage
import random
import math

class TestStorage(unittest.TestCase):

    def setUp(self):
        # Generate 1000 DNA sequences, of variying lenght and then create double stranded Oligos
        data = [''.join(random.choices(['A', 'C', 'G', 'T'], k=int(random.gauss(32*40, 20)))) for _ in range(1000)]
        data = [Oligo(sequence=ds) for ds in data]
        assembled_data = AssembledData()
        assembled_data.data = data
        return assembled_data


    def test_storage(self):
        # Simulate storage of DNA sequences in a storage medium
        data = self.setUp()

        storage_biogene = SimulateStorage(100, 'biogene')
        storage_permafrost = SimulateStorage(100, 'permafrost')

        data_stored_biogene = storage_biogene.simulate(data)
        data_stored_permafrost = storage_permafrost.simulate(data)   

        self.assertGreaterEqual(len(data.data), len(data_stored_biogene.data))
        self.assertIsInstance(data_stored_biogene, StoredData)
        self.assertGreaterEqual(len(data.data), len(data_stored_permafrost.data))
        self.assertIsInstance(data_stored_permafrost, StoredData)

if __name__ == '__main__':
    unittest.main()



