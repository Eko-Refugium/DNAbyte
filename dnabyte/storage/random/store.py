import random
import math
from dnabyte.store import SimulateStorage

def attributes(params):
    if 'years' not in params.__dict__ or params.years is None:
        years = 100
    else:
        years = params.years
        
    return {"years": years}

class Random(SimulateStorage):
    # Random: 1 cut/century/100 000 nucleotides
    # This class is for testing purposes mainly.
    def simulate(self, assembled_data):

        data_sto = []
        for oligo in assembled_data:
            for i in range(0, math.ceil(len(oligo[0])/460)):
                bases = ['A', 'T', 'C', 'G']
                randomposition = random.randint(0, len(oligo[0])-1)
                randombase = random.choice(bases)
                oligo[0] = oligo[0][:randomposition] + randombase + oligo[0][randomposition+1:]
            data_sto.append([oligo[0], oligo[1]])
            strand_breaks = math.ceil(len(oligo[0])/460)

        info = {"number_of_strand_breaks": strand_breaks}
        return data_sto, info
