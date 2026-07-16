import random
from dnabyte.sequence import SimulateSequencing

class IID(SimulateSequencing):
    """
    Simulate sequencing errors using an independent and identically distributed (IID) model.
    This model assumes that each base in the DNA sequence has a fixed probability of being
    substituted, inserted, or deleted, independent of the other bases and is for testing purposes only.
    """

    def simulate(self, data):
        """
        Simulate sequencing errors using an IID model.
        Supports substitutions, insertions, and deletions.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        # Get error parameters
        error_rate = getattr(self.params, 'iid_error_rate', 0.01)
        sub_rate = getattr(self.params, 'iid_substitution_rate', error_rate)
        ins_rate = getattr(self.params, 'iid_insertion_rate', 0.0)
        del_rate = getattr(self.params, 'iid_deletion_rate', 0.0)
        
        sequenceserror = []
        error_counter = 0
        
        for sequence in data:
            new_seq = list(sequence)
            i = 0
            
            while i < len(new_seq):
                rand_val = random.random()
                
                # Substitution error
                if rand_val < sub_rate:
                    error_counter += 1
                    new_base = random.choice(['A', 'C', 'G', 'T'])
                    new_seq[i] = new_base
                    i += 1
                    
                # Insertion error
                elif rand_val < sub_rate + ins_rate:
                    error_counter += 1
                    insert_base = random.choice(['A', 'C', 'G', 'T'])
                    new_seq.insert(i, insert_base)
                    i += 2  # Skip the inserted base
                    
                # Deletion error
                elif rand_val < sub_rate + ins_rate + del_rate:
                    error_counter += 1
                    new_seq.pop(i)
                    # Don't increment i, check next base at same position
                    
                else:
                    i += 1
            
            sequenceserror.append(''.join(new_seq))
               
        info = {
            'error_counter': error_counter,
            'substitution_rate': sub_rate,
            'insertion_rate': ins_rate,
            'deletion_rate': del_rate,
        }

        return sequenceserror, info
    
def check_parameter(parameter, default, min, max, inputparams):
    if not hasattr(inputparams, parameter) or inputparams.__dict__[parameter] is None:
        parameter_value = default
    elif not (min <= inputparams.__dict__[parameter] <= max):
        raise ValueError(f"{parameter} must be greater than or equal to {min} and less than or equal to {max}, got {inputparams.__dict__[parameter]}")
    else:
        parameter_value = inputparams.__dict__[parameter]
    
    return parameter_value
    
def attributes(params):
    iid_error_rate = check_parameter(
        parameter="iid_error_rate",
        default=0.01,
        min=0.0,
        max=1.0,
        inputparams=params
    )
    
    # Substitution rate (default to total error rate if not specified)
    iid_substitution_rate = check_parameter(
        parameter="iid_substitution_rate",
        default=iid_error_rate,
        min=0.0,
        max=1.0,
        inputparams=params
    )
    
    # Insertion rate (default 0 if not specified)
    iid_insertion_rate = check_parameter(
        parameter="iid_insertion_rate",
        default=0.0,
        min=0.0,
        max=1.0,
        inputparams=params
    )
    
    # Deletion rate (default 0 if not specified)
    iid_deletion_rate = check_parameter(
        parameter="iid_deletion_rate",
        default=0.0,
        min=0.0,
        max=1.0,
        inputparams=params
    )
        
    return {
        "iid_error_rate": iid_error_rate,
        "iid_substitution_rate": iid_substitution_rate,
        "iid_insertion_rate": iid_insertion_rate,
        "iid_deletion_rate": iid_deletion_rate,
    }