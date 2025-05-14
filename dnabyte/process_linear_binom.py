import random

from dnabyte.encode import Encode
from dnabyte.oligo import create_positional_libraries
from dnabyte.auxiliary import split_string_alternating_lengths, complementmap, transpose_matrix, sort_lists_by_first_n_entries, reduce_to_n_most_used_elements, find_closest_string

class ProcessLinearBinom(Encode):
    """
    This class implements data preprosessing for the binomial encoding before decoding.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    #def process_simple_binom(library, data.data, self.params.dna_barcode_length, **kwargs):
    def process_linear_binom(self, data):

        oligolength = len(self.params.library.library[0])
        positionallib = create_positional_libraries(self.params.library, len(self.params.library.leftmotives) + len(self.params.library.rightmotives) - 1)

        separated = []

        # TODO: parallelize
        for i in range(len(data.data)):
            fivetotreeprimetuples = []
            #for singular_codeword in codewords_tuples:
            fivetothreeend = data.data[i][0]
            
            lengthofcorrect = (len(self.params.library.leftmotives)+2*len(self.params.library.rightmotives) - 1) * oligolength 
            bases = ['A', 'C', 'T', 'G']
            
            while len(fivetothreeend) < lengthofcorrect:
                # Randomly choose a position to insert a base
                insert_position = random.randint(0, len(fivetothreeend))
                # Randomly choose a base to insert
                base_to_insert = random.choice(bases)
                # Insert the base at the chosen position
                fivetothreeend = fivetothreeend[:insert_position] + base_to_insert + fivetothreeend[insert_position:]
            
            while len(fivetothreeend) > lengthofcorrect:
                # Randomly choose a position to delete a base
                delete_position = random.randint(0, len(fivetothreeend) - 1)
                # Delete the base at the chosen position
                fivetothreeend = fivetothreeend[:delete_position] + fivetothreeend[delete_position + 1:]
                
            listofinfooligos = split_string_alternating_lengths(fivetothreeend, oligolength, oligolength // 2)
            
            for j in range(len(listofinfooligos)):
                if j % 2 == 1:
                    listofinfooligos[j] = (listofinfooligos[j][:len(listofinfooligos[j]) // 2][::-1] + listofinfooligos[j][len(listofinfooligos[j]) // 2:][::-1])[::-1]
                else:
                    listofinfooligos[j] = complementmap(listofinfooligos[j][:len(listofinfooligos[j])//2]) + complementmap(listofinfooligos[j][len(listofinfooligos[j]) // 2:])

            correctedlistofinfooligos = []
            
            for j, singularinformationoligo in enumerate(listofinfooligos):
                correctedlistofinfooligos.append(find_closest_string(singularinformationoligo, positionallib[j]))

            for j in range(self.params.dna_barcode_length):
                if correctedlistofinfooligos[j] in positionallib[j]:
                    correctedlistofinfooligos[j] = positionallib[j].index(correctedlistofinfooligos[j])
                # else:
                #     correctedlistofinfooligos[j]=random.randint(0,len(positionallib[j])-1)
            separated.append(correctedlistofinfooligos)
        listssorted = sort_lists_by_first_n_entries(separated, self.params.dna_barcode_length,theory = self.params.theory) 
        
        for i in range(len(listssorted)):
            for j in range(len(listssorted[i])):
                for k in range(self.params.dna_barcode_length):
                    listssorted[i][j].pop(0)
        
        listoflikley = []
        for evrycodeword in listssorted:
            listoflikley.append(reduce_to_n_most_used_elements(evrycodeword, self.params.sigma_amount))

        transposeinfo = []
        for i in range(len(listoflikley)):
            transposeinfo.append(transpose_matrix(listoflikley[i]))
        
        info = {}

        return transposeinfo, info