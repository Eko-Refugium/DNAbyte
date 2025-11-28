import numpy as np

from dnabyte.synthesize import SimulateSynthesis

class NoSynthPoly(SimulateSynthesis):


    def simulate(self, data):
        """
        Simulate synthesis of DNA sequences without any polymerase errors. This method is used for testing purposes.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        
        sequences_multiples = [seq for seq in data for _ in range(max(1, int(np.random.normal(self.mean, self.std_dev))))]
        average_copy_number = len(sequences_multiples) / len(data)
        info = {
            "average_copy_number": average_copy_number,
            "number_of_synthesis_errors": 0,
            "error_dict": {}
        }

        return sequences_multiples, info






