import random
import itertools

from dnabyte.encode import Encode
from dnabyte.auxiliary import dna_to_binary, bitstring_to_base4_list, sort_lists_by_first_n_entries_synth, most_common_string


class ProcessMaxDensity(Encode):
    """
    This class implements data preprosessing for the maxdensity encoding before decoding.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger
    
    bases = ['A', 'C', 'T', 'G']

    def process_maxdensity(self, data):

        # TODO: try to avoid spaces being in the data in the first place
        # Remove spaces from the data
        data.data = [seq.replace(' ', '') for seq in data.data] 
        
        # TODO: Why is this necessary? Can this be done in a better way?
        # Step 1: ensure that all codewords have the same length by filling random bases or deleting random bases
        for i in range(len(data.data)):
            

            while len(data.data[i]) < self.params.codeword_length:
                # Randomly choose a position to insert a base
                insert_position = random.randint(0, len(data.data[i]))
                # Randomly choose a base to insert
                base_to_insert = random.choice(self.bases)
                # Insert the base at the chosen position
                data.data[i] = data.data[i][:insert_position] + base_to_insert + data.data[i][insert_position:]
            
            while len(data.data[i]) > self.params.codeword_length:
                # Randomly choose a position to delete a base
                delete_position = random.randint(0, len(data.data[i]) - 1)
                # Delete the base at the chosen position
                data.data[i] = data.data[i][:delete_position] + data.data[i][delete_position + 1:]

        processed_list = []

        # TODO: parallelize
        for i in range(len(data.data)):
            
            # extract the barcode, which includes the index in the data object
            barcode = data.data[i][:self.params.dna_barcode_length]
            index_binary = dna_to_binary(barcode)
            index_base4_list = bitstring_to_base4_list(index_binary)
            index_base4_list.reverse()

            # TODO: Why is it necessary to save this as a list?
            indexed_list = []
            indexed_list.append(data.data[i][self.params.dna_barcode_length:])
            for j in range(len(index_base4_list)):
                indexed_list.insert(0, index_base4_list[j])

            lengthofthefirst = len(index_base4_list)
            processed_list.append(indexed_list)

        # Sort the list by the index numbers
        sorted_list = sort_lists_by_first_n_entries_synth(processed_list, lengthofthefirst) 

        # Remove the index numbers from the list
        for i in range(len(sorted_list)):
            for j in range(len(sorted_list[i])):
                for k in range(self.params.dna_barcode_length):
                    sorted_list[i][j].pop(0)
        
        # Step 2: find the most common string in each codeword
        list_of_most_common = []
        
        for codeword in sorted_list:
            list_of_most_common.append(most_common_string(codeword))
        
        info = {'number of codewords': len(list_of_most_common)}

        return list_of_most_common, info
        
