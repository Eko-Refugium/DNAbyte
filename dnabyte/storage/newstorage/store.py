from dnabyte.store import SimulateStorage
import random

class Newstorage(SimulateStorage):

    def __init__(self, params, logger=None):
        self.years = params.years

    def simulate(self, data):
        """
        Simulate storage of DNA sequences in cryogenic conditions.
        :param data: An object of class InSilicoDNA.
        :return: A tuple (stored_data, strand_breaks).
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in data:
            decay_probability = 1 - (1 - 1E-8)**(len(oligo) * self.years)
            if random.random() > decay_probability:
                remaining_oligos.append(oligo)
            else:
                strand_breaks += 1

        info = {"number_of_strand_breaks": strand_breaks}
        return remaining_oligos, info
    
def check_parameter(parameter, default, min, max, inputparams):
    if not hasattr(inputparams, parameter) or inputparams.__dict__[parameter] is None:
        parameter_value = default
    elif not (min <= inputparams.__dict__[parameter] <= max):
        raise ValueError(f"{parameter} must be greater than or equal to {min} and less than or equal to {max}, got {inputparams.__dict__[parameter]}")
    else:
        parameter_value = inputparams.__dict__[parameter]
    
    return parameter_value

def attributes(params):
    years = check_parameter(parameter="years",
                                        default=100,
                                        min=0,
                                        max=1000000000000,
                                        inputparams=params)

    return {"years": years}