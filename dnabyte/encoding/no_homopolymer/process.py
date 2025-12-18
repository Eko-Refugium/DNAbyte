import random
from collections import Counter
from tqdm import tqdm
from dnabyte.encoding.auxiliary import sort_lists_by_first_n_entries_synth

def process(data, params, logger=None):

    sortedlist = []

    # TODO: where are the spaces created?
    sequenceddata = [seq.replace(' ', '') for seq in data.data] 

    for i in tqdm(range(len(sequenceddata)), desc="Processing codewords", disable=logger is not None):
        
        bases = ['A', 'C', 'T', 'G']
        
        while len(sequenceddata[i]) < params.codeword_length:
            for tester2 in range(len(sequenceddata[i]) - 1):
                evenlist = ['T', 'G']
                oddlist = ['C', 'A']
                if sequenceddata[i][tester2] in evenlist and sequenceddata[i][tester2 + 1] in evenlist:
                    sequenceddata[i] = sequenceddata[i][:tester2 + 1] + random.choice(oddlist) + sequenceddata[i][tester2 + 2:]
                    break
                elif sequenceddata[i][tester2] in oddlist and sequenceddata[i][tester2 + 1] in oddlist:
                    sequenceddata[i] = sequenceddata[i][:tester2 + 1] + random.choice(evenlist) + sequenceddata[i][tester2 + 2:]
                    break
            else:
                # Randomly choose a position to insert a base
                insert_position = random.randint(0, len(sequenceddata[i]))
                # Randomly choose a base to insert
                base_to_insert = random.choice(bases)
                # Insert the base at the chosen position
                sequenceddata[i] = sequenceddata[i][:insert_position] + base_to_insert + sequenceddata[i][insert_position:]
                
        while len(sequenceddata[i]) > params.codeword_length:
            # Randomly choose a position to delete a base
            delete_position = random.randint(0, len(sequenceddata[i]) - 1)
            # Delete the base at the chosen position
            sequenceddata[i] = sequenceddata[i][:delete_position] + sequenceddata[i][delete_position + 1:]

        for tester2 in range(len(sequenceddata[i])):
            evenlist = ['T', 'G']
            oddlist = ['C', 'A']
            if tester2 % 2 == 0:
                if sequenceddata[i][tester2] not in evenlist:
                    sequenceddata[i] = sequenceddata[i][:tester2] + random.choice(evenlist) + sequenceddata[i][tester2 + 1:]
            else:
                if sequenceddata[i][tester2] not in oddlist:
                    sequenceddata[i] = sequenceddata[i][:tester2] + random.choice(oddlist) + sequenceddata[i][tester2 + 1:]

        indexcontainingdna = sequenceddata[i][:params.dna_barcode_length]
        indexbinary = dna_to_binary_custom(indexcontainingdna)
        index2list = []

        for j in range(len(indexbinary)):
            index2list.append(int(indexbinary[j]))
        index2list.reverse()
        sortablelist = []
        sortablelist.append(sequenceddata[i][params.dna_barcode_length:])
        for j in range(len(index2list)):
            sortablelist.insert(0, index2list[j])                    
        
        lengthofthefirst = len(index2list)
        sortedlist.append(sortablelist)
    
    listssorted = sort_lists_by_first_n_entries_synth(sortedlist, lengthofthefirst) 
    for i in range(len(listssorted)):
        for j in range(len(listssorted[i])):
            for k in range(params.dna_barcode_length):
                listssorted[i][j].pop(0)
    
    listoflikley = []
    
    for evrycodeword in tqdm(listssorted, desc="Finding consensus", disable=logger is not None):
        listoflikley.append(most_common_string(evrycodeword))

    info = {}

    return listoflikley, info

def dna_to_binary_custom(dna_string):
    binary_string = []
    for i, base in enumerate(dna_string):
        if i % 2 == 0:  # Even position
            if base == 'T':
                binary_string.append('0')
            elif base == 'G':
                binary_string.append('1')
        else:  # Odd position
            if base == 'C':
                binary_string.append('0')
            elif base == 'A':
                binary_string.append('1')
    return ''.join(binary_string)

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

