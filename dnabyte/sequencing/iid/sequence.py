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
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        sequenceserror = list(data)  # Create a copy of the data
        error_counter = 0
        for i in range(len(sequenceserror)):
            for k in range(len(sequenceserror[i])):
                if random.random() < self.params.iid_error_rate:
                    error_counter += 1
                    # randomly select a new base
                    new_base = random.choice(['A', 'C', 'G', 'T'])
                    # replace the base at the mutation position with the new base
                    sequenceserror[i] = (sequenceserror[i][:k] + new_base + sequenceserror[i][k + 1:])      
               
        info = {}
        info['error_counter'] = error_counter

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
    iid_error_rate = check_parameter(parameter="iid_error_rate",
                                        default=0.01,
                                        min=0.0,
                                        max=1.0,
                                        inputparams=params)
        
    return {"iid_error_rate": iid_error_rate}