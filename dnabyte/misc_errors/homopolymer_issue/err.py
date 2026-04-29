import random
from dnabyte.misc_err import SimulateMiscErrors

class HomopolymerIssue(SimulateMiscErrors):
    """
    Simulate sequencing errors using an independent and identically distributed (IID) model.
    This model assumes that each base in the DNA sequence has a fixed probability of being
    substituted, inserted, or deleted, independent of the other bases and is for testing purposes only.
    """

    def __init__(self, params):
        if not hasattr(params, 'deletion_prob') or params.deletion_prob is None:
            self.deletion_prob = 0.01
        else:
            self.deletion_prob = params.deletion_prob

    def simulate(self, data):
        """
        Simulate sequencing errors using a homopolymer issue model.
        
        :param data: A list of DNA sequences (strings).
        :return: A list of sequenced DNA sequences.
        """
        triplets = ['AAAA', 'CCCC', 'TTTT', 'GGGG']
        sequenceserror = []
        error_counter = 0
        for seq in data:
            sequence = list(seq)
            ticker = 0
            while ticker <= len(sequence) - 4:
                triplet = ''.join(sequence[ticker:ticker+4])
                if triplet in triplets:
                    if random.random() < self.deletion_prob:
                        # Delete one random letter from the triplet
                        del_index = ticker + random.randint(0, 2)
                        sequence.pop(del_index)
                        ticker = max(ticker - 2, 0)
                        error_counter += 1
                        continue  # Skip increment to handle updated indices
                ticker += 1
            sequenceserror.append(''.join(sequence))
        info = {'error_counter': error_counter}
        return sequenceserror, info
    
def attributes(params):
    if 'deletion_prob' not in params.__dict__ or params.deletion_prob is None:
        deletion_prob = 0.01
    else:
        deletion_prob = params.deletion_prob
        
    return {"deletion_prob": deletion_prob}