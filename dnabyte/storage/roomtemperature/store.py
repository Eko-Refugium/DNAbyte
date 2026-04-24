import random
from dnabyte.store import SimulateStorage

class Roomtemperature(SimulateStorage):
    # Reference::
    # half-life of DNA: 521 years => lambda = 0.00133
    
    def __init__(self, params, logger=None):
        self.years = params.years

    def simulate(self, data):
        """
        Simulate storage of DNA sequences at room temperature.
        :param assembled_data: An object of class AssembledData.
        :return: A list of stored DNA sequences.
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in data:
            # half-life of DNA: 521 years => lambda = ln(2)/521 ≈ 0.00133 per year
            # decay rate per nucleotide per year
            lambda_decay = 0.00133
            decay_probability = 1 - (1 - lambda_decay)**(len(oligo) * self.years)
            if random.random() > decay_probability:
                remaining_oligos.append(oligo)
            else:
                strand_breaks += 1

        info = {"number_of_strand_breaks": strand_breaks}
        if remaining_oligos == []:
            raise ValueError("All DNA strands have decayed during storage at room temperature. No sequences remain.")
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