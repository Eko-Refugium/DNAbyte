import random

from dnabyte.encode import Encode
from dnabyte.auxiliary import split_string_into_chunks_poly_chain, find_closest_string, sort_lists_by_first_n_entries, reduce_to_n_most_used_elements

class ProcessPolyChain(Encode):
    """
    This class implements data preprosessing for the linear encoding before decoding.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger
            
    def process_poly_chain(self, data):

        DNAs = sorted(self.params.library.messages)

        position = self.params.library.position
        position_length = len(position[0])
        oligo_length = len(DNAs[0])
        generic_length = len(self.params.library.generic[0])
        
        separated = []

        for i in range(len(data.data)):
            
            fivetothreeend = data.data[i][0]

            if len(position) % 2 == 0:
                lengthofcorrect = (len(position) - 1) * (oligo_length + 2 * generic_length) + (len(position)) * position_length
            else:
                lengthofcorrect = (len(position) - 1) * (oligo_length + 2 * generic_length) + (len(position) - 1) * position_length
            bases = ['A', 'C', 'T', 'G']

            while len(fivetothreeend) < lengthofcorrect:
                # Randomly choose a position to insert a base
                insert_position = random.randint(self.params.dna_barcode_length, len(fivetothreeend))
                # Randomly choose a base to insert
                base_to_insert = random.choice(bases)
                # Insert the base at the chosen position
                fivetothreeend = fivetothreeend[:insert_position] + base_to_insert + fivetothreeend[insert_position:]
            
            while len(fivetothreeend) > lengthofcorrect:
                # Randomly choose a position to delete a base
                delete_position = random.randint(0, len(fivetothreeend) - 1)
                # Delete the base at the chosen position
                fivetothreeend = fivetothreeend[:delete_position] + fivetothreeend[delete_position + 1:]
            
            chunks = split_string_into_chunks_poly_chain(fivetothreeend,generic_length,position_length,oligo_length)
            correctedlistofinfooligos = []
            for singularinformationoligo in chunks:
                correctedlistofinfooligos.append(find_closest_string(singularinformationoligo, DNAs))
            for j in range(self.params.dna_barcode_length):
                correctedlistofinfooligos[j] = DNAs.index(correctedlistofinfooligos[j])
            separated.append(correctedlistofinfooligos)
        
        listssorted = sort_lists_by_first_n_entries(separated, self.params.dna_barcode_length, theory = self.params.theory) 

        for i in range(len(listssorted)):
            for j in range(len(listssorted[i])):
                for k in range(self.params.dna_barcode_length):
                    listssorted[i][j].pop(0)

        listoflikley = []
        for evrycodeword in listssorted:
            listoflikley.append(reduce_to_n_most_used_elements(evrycodeword,1))
        for i in range(len(listoflikley)):
            listoflikley[i] = listoflikley[i][0]
        
        info = {}

        return listoflikley, info