import unittest
from dnabyte.store import SimulateStorage
from backup.data import AssembledData
import random

class TestSimulateStorage(unittest.TestCase):
    storage_conditions = ['permafrost', 'room_temperature', 'biogene']

    def setUp(self):
        # Set parameters
        self.k = 32 * 40  # Length of oligo to be sequenced
        self.n = 1000  # Number of oligos to be sequenced
        
        # Create random DNA sequences
        self.data = [''.join(random.choices(['A', 'C', 'G', 'T'], k=self.k)) for _ in range(self.n)]
        # Create StoredData object
        self.data_before = AssembledData(encoded_data=self.data)

    def test(self):

        for condition in self.storage_conditions:

            # logger.info('Storage condition: ' + condition)
            # Instantiate SimulateSequencing object
            simulator = SimulateStorage(storage_condition=condition)

            # Simulate sequencing
            data_after = simulator.simulate(self.data_before)

            # Count differences
            cntr = 0
            for i in range(len(self.data_before.data)):
                if self.data_before.data[i] != data_after.data[i]:
                    cntr += 1

            # Calculate percentage of differences
            percentage_of_differences = cntr / len(self.data_before.data) * 100
            logger.info('Percentage of differences:', percentage_of_differences)

            # Assert that there are differences
            self.assertGreater(percentage_of_differences, -1, f"Storage Condition: {condition}")


if __name__ == '__main__':
    unittest.main()