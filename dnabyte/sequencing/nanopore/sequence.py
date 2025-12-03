import numpy as np

from dnabyte.sequence import SimulateSequencing

class Nanopore(SimulateSequencing):
    """
    Nanopore sequencing errors are characterized by a higher error rate compared to Illumina sequencing, with a mix of substitution, 
    insertion, and deletion errors. A good model for simulating substitution errors in Nanopore sequencing should account for these characteristics:
    Higher Error Rates: Nanopore sequencing typically has higher error rates, often around 5-15%.
    Error Distribution: Errors are more uniformly distributed across the read length compared to Illumina sequencing.
    Substitution Matrix: Use a substitution matrix to define the probabilities of substituting one base for another.
    """

    def simulate(self, data):
        """
        Simulate sequencing errors using Illumina technology.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """

        substitution_matrix_nanopore = {
                'A': {'A': 0.85, 'T': 0.05, 'C': 0.05, 'G': 0.05},
                'T': {'A': 0.05, 'T': 0.85, 'C': 0.05, 'G': 0.05},
                'C': {'A': 0.05, 'T': 0.05, 'C': 0.85, 'G': 0.05},
                'G': {'A': 0.05, 'T': 0.05, 'C': 0.05, 'G': 0.85}
            }

        sequenced_data = []
        error_counter = 0

        for sequence in data.data:

            # Simulate substitutions
            modified_sequence, counter = self.simulate_substitutions(sequence, substitution_matrix_nanopore, k=0.01)
            error_counter += counter

            # Simulate instertions
            modified_sequence, counter = self.simulate_insertions(modified_sequence, ins_prob=0.01)
            error_counter += counter

            # Simulate deletions
            modified_sequence, counter = self.simulate_deletions(modified_sequence, del_prob=0.01)
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


def attributes(params):
    return {}


