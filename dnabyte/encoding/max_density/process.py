import random
from collections import Counter

from dnabyte.encoding.auxiliary import sort_lists_by_first_n_entries_synth

def process(data, params, logger=None):

    bases = ['A', 'C', 'T', 'G']

    # TODO: try to avoid spaces being in the data in the first place
    # Remove spaces from the data
    data.data = [seq.replace(' ', '') for seq in data.data] 
    
    # TODO: Why is this necessary? Can this be done in a better way?
    # Step 1: ensure that all codewords have the same length by filling random bases or deleting random bases
    for i in range(len(data.data)):
        

        while len(data.data[i]) < params.codeword_length:
            # Randomly choose a position to insert a base
            insert_position = random.randint(0, len(data.data[i]))
            # Randomly choose a base to insert
            base_to_insert = random.choice(bases)
            # Insert the base at the chosen position
            data.data[i] = data.data[i][:insert_position] + base_to_insert + data.data[i][insert_position:]
        
        while len(data.data[i]) > params.codeword_length:
            # Randomly choose a position to delete a base
            delete_position = random.randint(0, len(data.data[i]) - 1)
            # Delete the base at the chosen position
            data.data[i] = data.data[i][:delete_position] + data.data[i][delete_position + 1:]

    processed_list = []

    # TODO: parallelize
    for i in range(len(data.data)):
        
        # extract the barcode, which includes the index in the data object
        barcode = data.data[i][:params.dna_barcode_length]
        index_binary = dna_to_binary(barcode)
        index_base4_list = bitstring_to_base4_list(index_binary)
        index_base4_list.reverse()

        # TODO: Why is it necessary to save this as a list?
        indexed_list = []
        indexed_list.append(data.data[i][params.dna_barcode_length:])
        for j in range(len(index_base4_list)):
            indexed_list.insert(0, index_base4_list[j])

        lengthofthefirst = len(index_base4_list)
        processed_list.append(indexed_list)

    # Sort the list by the index numbers
    sorted_list = sort_lists_by_first_n_entries_synth(processed_list, lengthofthefirst) 

    # Remove the index numbers from the list
    for i in range(len(sorted_list)):
        for j in range(len(sorted_list[i])):
            for k in range(params.dna_barcode_length):
                sorted_list[i][j].pop(0)
    
    # Step 2: find the most common string in each codeword
    list_of_most_common = []
    
    for codeword in sorted_list:
        list_of_most_common.append(most_common_string(codeword))
    
    info = {'number of codewords': len(list_of_most_common)}

    return list_of_most_common, info
    
# TODO: We need to implement a few sanity checks here
def most_common_string(strings):
    if not strings:
        return ""
    
    # Initialize a list to store the most common letters at each position
    most_common_letters = []
    
    # Iterate over each position in the strings
    for i in range(len(strings[0])):
        # Get the letters at the current position from all strings
        letters_at_position = [s[i] for s in strings]
        
        # Find the most common letter at the current position
        most_common_letter = Counter(letters_at_position).most_common(1)[0][0]
        
        # Append the most common letter to the list
        most_common_letters.append(most_common_letter)
    
    # Join the list of most common letters into a single string
    return ''.join(most_common_letters)

def bitstring_to_base4_list(bitstring):
    base4_list = [int(bitstring[i:i+2], 2) for i in range(0, len(bitstring), 2)]
    return base4_list

def dna_to_binary(dna_string):
    
    binary_mapping = {
        'A': '00',
        'G': '01',
        'C': '10',
        'T': '11'
    }
    binary_string = ''.join(binary_mapping[base] for base in dna_string)
    return binary_string
