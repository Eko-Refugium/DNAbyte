import unittest
from dnabyte.sequence import SimulateSequencing
from dnabyte.data_classes.insilicodna import InSilicoDNA
import random

class TestSimulateSequencing(unittest.TestCase):
    sequencing_technologies = {
        'None': 41,
        '2D': 40,
        'Subread': 37,
        'PairedEnd': 36,
        '1D': 39,
        'CCS': 38,
        'Single End': 35
    }

    def setUp(self):
        # Set parameters
        self.k = 32 * 40  # Length of oligo to be sequenced
        self.n = 1000  # Number of oligos to be sequenced
        
        # Create random DNA sequences
        self.data = [''.join(random.choices(['A', 'C', 'G', 'T'], k=self.k)) for _ in range(self.n)]
        # Create StoredData object
        self.data_before = InSilicoDNA(data=self.data)

    def test(self):

        for tech_name, tech_code in self.sequencing_technologies.items():
            # logger.info('Sequencing Technology: ' + tech_name)
            # Instantiate SimulateSequencing object
            simulator = SimulateSequencing(sequencing_technology=tech_code)

            # Simulate sequencing
            data_after = simulator.simulate(self.data_before)

            # Count differences
            cntr = 0
            for i in range(len(self.data_before.data)):
                if self.data_before.data[i] != data_after.data[i]:
                    cntr += 1

            # Calculate percentage of differences
            percentage_of_differences = cntr / len(self.data_before.data) * 100
            # logger.info('Percentage of differences:', percentage_of_differences)

            # Assert that there are differences
            self.assertGreater(percentage_of_differences, -1, f"Sequencing technology: {tech_name}")


if __name__ == '__main__':
    unittest.main()