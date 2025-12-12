import random
import math
from dnabyte.store import SimulateStorage


class Random(SimulateStorage):
    # Random: 1 cut/century/100 000 nucleotides
    # This class is for testing purposes mainly.

    def __init__(self, params, logger=None):
        self.years = params.years

    def simulate(self, data):

        remaining_oligos = []
        total_strand_breaks = 0
        
        for oligo in data:
            for i in range(0, math.ceil(len(oligo[0])/460)):
                bases = ['A', 'T', 'C', 'G']
                randomposition = random.randint(0, len(oligo[0])-1)
                randombase = random.choice(bases)
                oligo[0] = oligo[0][:randomposition] + randombase + oligo[0][randomposition+1:]
            remaining_oligos.append([oligo[0], oligo[1]])
            strand_breaks = math.ceil(len(oligo[0])/460)

        info = {"number_of_strand_breaks": total_strand_breaks}
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