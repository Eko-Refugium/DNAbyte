import numpy as np
from dnabyte.data import SequencedData, StoredData
from dnabyte.ErrorChannels.sequencing_simulation import sequencing_simulation

class SimulateSequencing:
    """
    Simulate sequencing of DNA sequences.
    There are two main types of sequencing technologies: Illumina and Nanopore, each with their 
    characteristic error profiles.

    :param sequencing_technology: A string 'illumina', or 'nanopore'. None to perform no sequencing simulation.
    """
    
    def __init__(self, params, logger=None):
        self.sequencing_method = params.sequencing_method
        if self.sequencing_method == 'iid':
            self.iidmethod = params.iid_error_rate
        else:
            self.iidmethod = None
        
    def simulate(self, data):
        """
        Simulate sequencing of DNA sequences.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        if self.sequencing_method == None:
            sequenced_data = SequencedData(data=data.data)
            error_dict = {}
            
        else:
            if isinstance(data, StoredData):
                
                # # Simulate sequencing errors
                # if self.sequencing_technology == 'illumina':
                #     sequenced_data = [self.simulate_illumina_sequencing(seq) for seq in data]
                # elif self.sequencing_technology == 'nanopore':
                #     sequenced_data = [self.simulate_nanopore_sequencing(seq) for seq in data]
                # else:
                #     raise ValueError("Sequencing technology not recognized.")
                
                data_seq, error_dict = sequencing_simulation(data.data, method_id=self.sequencing_method, errorrate=self.iidmethod)
                
                sequenced_data = SequencedData(data=data_seq)
            else:
                raise ValueError("The input data is not an instance of StoredData.")
        info = {}
        info['eror_dict'] = error_dict

        return sequenced_data, info


    # # Define a substitution matrix
    # substitution_matrix_illumina = {
    #     'A': {'A': 0.997, 'T': 0.001, 'C': 0.0005, 'G': 0.0015},
    #     'T': {'A': 0.001, 'T': 0.997, 'C': 0.0015, 'G': 0.0005},
    #     'C': {'A': 0.0005, 'T': 0.0015, 'C': 0.007, 'G': 0.001},
    #     'G': {'A': 0.0015, 'T': 0.0005, 'C': 0.001, 'G': 0.997}
    # }

    # substitution_matrix_nanopore = {
    #         'A': {'A': 0.85, 'T': 0.05, 'C': 0.05, 'G': 0.05},
    #         'T': {'A': 0.05, 'T': 0.85, 'C': 0.05, 'G': 0.05},
    #         'C': {'A': 0.05, 'T': 0.05, 'C': 0.85, 'G': 0.05},
    #         'G': {'A': 0.05, 'T': 0.05, 'C': 0.05, 'G': 0.85}
    #     }

    # def simulate_illumina_sequencing(self, sequence):
    #     """
    #     Illumina errors are not random, but are strand specific and the error rates are higher at the end of a read.
    #     Substitution error rates are about 0.0015 - 0.0004 errors per base.
    #     Insertions and deletions are significantly less likely, they are on the order of 10^-6.
    #     """

    #     modified_sequence = self.simulate_substitutions(sequence, self.substitution_matrix_illumina, k=0.01)
    #     modified_sequence = self.simulate_insertions(modified_sequence, ins_prob=1e-6)
    #     modified_sequence = self.simulate_deletions(modified_sequence, del_prob=1e-6)

    #     return modified_sequence


    # def simulate_nanopore_sequencing(self, sequence):
    #     """
    #     Nanopore sequencing errors are characterized by a higher error rate compared to Illumina sequencing, with a mix of substitution, 
    #     insertion, and deletion errors. A good model for simulating substitution errors in Nanopore sequencing should account for these characteristics:

    #     Higher Error Rates: Nanopore sequencing typically has higher error rates, often around 5-15%.
    #     Error Distribution: Errors are more uniformly distributed across the read length compared to Illumina sequencing.
    #     Substitution Matrix: Use a substitution matrix to define the probabilities of substituting one base for another.
    #     """

    #     modified_sequence = self.simulate_substitutions(sequence, self.substitution_matrix_illumina_nanopore, k=0.01)
    #     modified_sequence = self.simulate_insertions(modified_sequence, ins_prob=0.01)
    #     modified_sequence = self.simulate_deletions(modified_sequence, del_prob=0.01)

    #     return modified_sequence
    

    # # Simulate substitutions in a sequence
    # def simulate_substitutions(self, sequence, substitution_matrix_illumina, k=0.01):
    #     L = len(sequence)
    #     new_sequence = list(sequence)
    #     for p in range(L):
    #         if np.random.rand() < self.positional_error_rate(p, L, k):
    #             original_base = sequence[p]
    #             new_sequence[p] = self.get_new_base(original_base, substitution_matrix_illumina)
    #     return ''.join(new_sequence)

    # # Positional error rate function
    # def positional_error_rate(p, L, k=0.01):
    #     """
    #     This function enables to model the error rate as a function of the position in the read.
    #     """    
    #     return k * (p / L)


    # def simulate_insertions(self, sequence, ins_prob=1e-6):
    #     """
    #     Simulate insertions in a sequence.
    #     """
    #     L = len(sequence)
    #     new_sequence = list(sequence)
    #     for p in range(L):
    #         if np.random.rand() < ins_prob:
    #             new_sequence.insert(p, np.random.choice(['A', 'T', 'C', 'G']))
    #     return ''.join(new_sequence)
    

    # def simulate_deletions(self, sequence, del_prob=1e-6):
    #     """
    #     Simulate deletions in a sequence.
    #     """
    #     L = len(sequence)
    #     new_sequence = list(sequence)
    #     for p in range(L):
    #         if np.random.rand() < del_prob:
    #             new_sequence.pop(p)
    #     return ''.join(new_sequence)
    

    # # Function to get a new base based on substitution probabilities
    # def get_new_base(self, original_base, substitution_matrix_illumina):
    #     bases = ['A', 'T', 'C', 'G']
    #     probabilities = [substitution_matrix_illumina[original_base][b] for b in bases]
    #     total_prob = sum(probabilities)
    #     probabilities = [p / total_prob for p in probabilities]  # Normalize to ensure sum is 1
    #     return np.random.choice(bases, p=probabilities)










