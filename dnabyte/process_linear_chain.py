import random

from dnabyte.encode import Encode
from dnabyte.auxiliary import split_string_into_chunks, find_closest_string, sort_lists_by_first_n_entries, reduce_to_n_most_used_elements


class ProcessLinearChain(Encode):
    """
    This class implements data preprosessing for the encoding before decoding.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger
        self.bases = ['A', 'C', 'T', 'G']
    def process_linear_chain(self, data):

        oligolength = len(self.params.library.library[0])
        separated = []
        
        # TODO: parallelize
        for i in range(len(data.data)):            

            if self.params.codeword_length % 2 == 0:
                lengthofcorrect = (self.params.codeword_length) * (len(self.params.library.library[0]))
            else:
                lengthofcorrect = self.params.codeword_length*(len(self.params.library.library[0]))
            
            while len(data.data[i][1]) < lengthofcorrect:
                # Randomly choose a position to insert a base
                insert_position = random.randint(0, len(data.data[i][1]))
                # Randomly choose a base to insert
                base_to_insert = random.choice(self.bases)
                # Insert the base at the chosen position
                data.data[i][1] = data.data[i][1][:insert_position] + base_to_insert + data.data[i][1][insert_position:]
            
            while len(data.data[i][1]) > lengthofcorrect:
                # Randomly choose a position to delete a base
                delete_position = random.randint(0, len(data.data[i][1]) - 1)
                # Delete the base at the chosen position
                data.data[i][1] = data.data[i][1][:delete_position] + data.data[i][1][delete_position + 1:]
                
            chunks = split_string_into_chunks(data.data[i][1],oligolength)
            correctedlistofinfooligos = []
            for singularinformationoligo in chunks:
                correctedlistofinfooligos.append(find_closest_string(singularinformationoligo, self.params.library.library))
            for j in range(self.params.dna_barcode_length):
                correctedlistofinfooligos[j]=self.params.library.library.index(correctedlistofinfooligos[j])
            separated.append(correctedlistofinfooligos)
            
        listssorted = sort_lists_by_first_n_entries(separated, self.params.dna_barcode_length, self.params.theory) 
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