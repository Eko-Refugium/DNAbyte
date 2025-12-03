import random
from dnabyte.store import SimulateStorage

class Permafrost(SimulateStorage):
    # permafrost: 5.5E−6/nt/yr
    # Reference:
    # Allentoft, M. E. et al. The half-life of DNA in bone: measuring decay kinetics in 158 dated fossils. 
    # Proc. R. Soc. B Biol. Sci. 279, 4724–4733 (2012).
    def simulate(self, assembled_data):
        """
        Simulate storage of DNA sequences in a permafrost medium.
        :param assembled_data: An object of class AssembledData.
        :return: A list of stored DNA sequences.
        """
        remaining_oligos = []
        strand_breaks = 0

        for oligo in assembled_data:
            decay_probability = 1 - (1 - 5.5E-6)**(len(oligo) * self.years)
            if random.random() > decay_probability:
                remaining_oligos.append(oligo)
            else:
                strand_breaks += 1

        info = {"number_of_strand_breaks": strand_breaks}
        return remaining_oligos, info
    
def attributes(params):
    if 'years' not in params.__dict__ or params.years is None:
        years = 100
    else:
        years = params.years
        
    return {"years": years}
