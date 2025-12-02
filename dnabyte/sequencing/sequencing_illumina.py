import numpy as np

from dnabyte.sequence import SimulateSequencing

def attributes(params):
    return {}

class Illumina(SimulateSequencing):
    """
    Simulate sequencing errors using Illumina technology.
    Illumina errors are not random, but are strand specific and the error rates are higher at the end of a read.
    Substitution error rates are about 0.0015 - 0.0004 errors per base.
    Insertions and deletions are significantly less likely, they are on the order of 10^-6.
    """

    def simulate(self, data):
        """
        Simulate sequencing errors using Illumina technology.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """

        substitution_matrix_illumina = {
        'A': {'A': 0.997, 'T': 0.001, 'C': 0.0005, 'G': 0.0015},
        'T': {'A': 0.001, 'T': 0.997, 'C': 0.0015, 'G': 0.0005},
        'C': {'A': 0.0005, 'T': 0.0015, 'C': 0.007, 'G': 0.001},
        'G': {'A': 0.0015, 'T': 0.0005, 'C': 0.001, 'G': 0.997}
        }  

        sequenced_data = []
        error_counter = 0

        for sequence in data.data:
            # Simulate substitutions
            modified_sequence, counter = self.simulate_substitutions(sequence, substitution_matrix_illumina, k=0.01)
            error_counter += counter

            # Simulate instertions
            modified_sequence, counter = self.simulate_insertions(modified_sequence, ins_prob=1e-6)
            error_counter += counter

            # Simulate deletions
            modified_sequence, counter = self.simulate_deletions(modified_sequence, del_prob=1e-6)
            error_counter += counter

            sequenced_data.append(modified_sequence)

        return sequenced_data, error_counter

    def simulate_substitutions(self, sequence, substitution_matrix_illumina, k=0.01):
        new_sequence = list(sequence)
        counter = 0
        for p in range(len(sequence)):
            if np.random.rand() < self.positional_error_rate(p, len(sequence), k):
                counter += 1
                original_base = sequence[p]
                new_sequence[p] = self.get_new_base(original_base, substitution_matrix_illumina)
        return ''.join(new_sequence), counter
    

    def simulate_insertions(self, sequence, ins_prob=1e-6):
        """
        Simulate insertions in a sequence.
        """
        new_sequence = list(sequence)
        counter = 0
        for p in range(len(sequence)):
            if np.random.rand() < ins_prob:
                counter += 1
                new_sequence.insert(p, np.random.choice(['A', 'T', 'C', 'G']))
        return ''.join(new_sequence), counter
    

    def simulate_deletions(self, sequence, del_prob=1e-6):
        """
        Simulate deletions in a sequence.
        """
        new_sequence = list(sequence)
        counter = 0
        for p in range(len(sequence)):
            if np.random.rand() < del_prob:
                counter += 1
                new_sequence.pop(p)
        return ''.join(new_sequence), counter
    
    # Positional error rate function
    def positional_error_rate(p, L, k=0.01):
        """
        This function enables to model the error rate as a function of the position in the read.
        """    
        return k * (p / L)

    # Function to get a new base based on substitution probabilities
    def get_new_base(self, original_base, substitution_matrix_illumina):
        bases = ['A', 'T', 'C', 'G']
        probabilities = [substitution_matrix_illumina[original_base][b] for b in bases]
        total_prob = sum(probabilities)
        probabilities = [p / total_prob for p in probabilities]  # Normalize to ensure sum is 1
        return np.random.choice(bases, p=probabilities)




