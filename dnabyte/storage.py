from dnabyte.data import StoredData, AssembledData
import random
import math

# room temperature:
# half-life of DNA: 521 years => lambda = 0.00133

# permafrost: 5.5E−6/nt/yr
# Reference:
# Allentoft, M. E. et al. The half-life of DNA in bone: measuring decay kinetics in 158 dated fossils. 
# Proc. R. Soc. B Biol. Sci. 279, 4724–4733 (2012).

# biogene: (anhydrous and anoxic atmosphere maintained inside hermetic capsules)
# 1E-7/nt/yr =?= 1 cut/century/100 000 nucleotides = 
# Reference: 
# Coudy, Delphine, et al. "Long term conservation of DNA at ambient temperature. 
# Implications for DNA data storage." PLoS One 16.11 (2021): e0259868.


class SimulateStorage:
    """
    Simulate storage of DNA sequences in a storage medium.

    The class consists of a constructor, which sets the parameters of the storage simulation, 
    and a simulate method, which applies the simulation to an object of the AssembledData class.
    
    :param years: The number of years for which DNA storage is supposed to be simulated. (0 to perform no simulation)
    :param storage_conditions: A string 'permafrost', or 'room_temperature', or 'biogene'. (None to perform no simulation)
    """
    def __init__(self, params, logger=None):

        self.years = params.years
        self.storage_conditions = params.storage_conditions

    def simulate(self, assembled_data):
        """
        Simulate storage of DNA sequences in a storage medium.
        
        :param assembled_data: An object of class AssembledData.
        :return: A list of stored DNA sequences.
        """
        
        if self.storage_conditions == None or self.years == 0:
            stored_data = StoredData(assembled_data=assembled_data.data)
            strand_breaks = 0

        elif self.storage_conditions == 'random':
            holder = []
            for oligo in assembled_data.data:
                for i in range(0, math.ceil(len(oligo[0])/460)):
                    bases = ['A', 'T', 'C', 'G']
                    randomposition = random.randint(0, len(oligo[0])-1)
                    randombase = random.choice(bases)
                    oligo[0] = oligo[0][:randomposition] + randombase + oligo[0][randomposition+1:]
                holder.append([oligo[0], oligo[1]])
                strand_breaks = math.ceil(len(oligo[0])/460)

            stored_data = StoredData(assembled_data=holder)
                
        else:
            if isinstance(assembled_data, AssembledData):
                data = assembled_data.data
                #scaled_data = [item for item in data for _ in range(self.scale)]
                remaining_oligos = []
                #for oligo in scaled_data:

                counter = 0

                for oligo in data:

                # Probability of decay per oligo equals 1 minus the probability that not a single nucleotide decayed
                    if self.storage_conditions == 'permafrost':
                        decay_probability = 1 - (1 - 5.5E-6)**(len(oligo)*self.years)
                    elif self.storage_conditions == 'biogene':
                        decay_probability = 1 - (1 - 1E-7)**(len(oligo)*self.years)

                    else: 
                        decay_probability = 1 - (1 - 5E-3)**(len(oligo)*self.years)
                                            
                    rndm = random.random()

                    if rndm > decay_probability:
                        remaining_oligos.append(oligo)
                    else:
                        counter += 1

                    strand_breaks = counter
                    # Approach, where every oligo has an equal probability of decay
                    # N_t = len(data) * math.exp(-decay_probability * self.years)
                    # remaining_oligos = random.sample(data, int(round(N_t)))
                                                    
                stored_data = StoredData(assembled_data=remaining_oligos)
                                                    
            else: 
                raise ValueError("The input data is not an instance of StoredData.")

        info = {}
        info['number_of_strand_breaks'] = strand_breaks

        return stored_data, info

