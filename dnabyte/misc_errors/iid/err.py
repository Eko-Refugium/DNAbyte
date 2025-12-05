import random
from dnabyte.misc_err import SimulateMiscErrors

class Err_IID(SimulateMiscErrors):
    """
    Simulate sequencing errors using an independent and identically distributed (IID) model.
    This model assumes that each base in the DNA sequence has a fixed probability of being
    substituted, inserted, or deleted, independent of the other bases and is for testing purposes only.
    """
    def __init__(self, params):
        if not hasattr(params, 'iid_error_rate') or params.iid_error_rate is None:
            self.iid_error_rate = 0.01
        else:
            self.iid_error_rate = params.iid_error_rate
            
    def simulate(self, data):
        """
        Simulate sequencing errors using an IID model.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        print('self if err:', self)
        sequenceserror = list(data)  # Create a copy of the data
        error_counter = 0
        for i in range(len(sequenceserror)):
            for j in range(len(sequenceserror[i])):
                for k in range(len(sequenceserror[i][j])):
                    if random.random() < self.iid_error_rate:
                        error_counter += 1
                        # randomly select a new base
                        new_base = random.choice(['A', 'C', 'G', 'T'])
                        # replace the base at the mutation position with the new base
                        sequenceserror[i][j] = (data[i][j][:k] + new_base + data[i][j][k + 1:])      
               
        info = {}
        info['error_counter'] = error_counter

        return sequenceserror, info
    
def attributes(params):
    if 'iid_error_rate' not in params.__dict__ or params.iid_error_rate is None:
        iid_error_rate = 0.01
    else:
        iid_error_rate = params.iid_error_rate
        
    return {"iid_error_rate": iid_error_rate}