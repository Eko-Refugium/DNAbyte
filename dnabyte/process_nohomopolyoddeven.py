import random

from dnabyte.encode import Encode
from dnabyte.auxiliary import dna_to_binary_custom, sort_lists_by_first_n_entries_synth, most_common_string


class ProcessNoHomoPolyOddEven(Encode):
    """
    This class implements data preprosessing for the nohomopolymer encoding before decoding.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def process_nohomopolyoddeven(self, data):

        sortedlist = []

        # TODO: where are the spaces created?
        sequenceddata = [seq.replace(' ', '') for seq in data.data] 

        for i in range(len(sequenceddata)):
            
            bases = ['A', 'C', 'T', 'G']
            
            while len(sequenceddata[i]) < self.params.codeword_length:
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
                    
            while len(sequenceddata[i]) > self.params.codeword_length:
                # Randomly choose a position to delete a base
                delete_position = random.randint(0, len(sequenceddata[i]) - 1)
                # Delete the base at the chosen position
                sequenceddata[i] = sequenceddata[i][:delete_position] + sequenceddata[i][delete_position + 1:]

            for tester2 in range(len(sequenceddata[i]) - 1):
                evenlist = ['T', 'G']
                oddlist = ['C', 'A']
                if tester2 % 2 == 0:
                    if sequenceddata[i][tester2] not in evenlist:
                        sequenceddata[i] = sequenceddata[i][:tester2] + random.choice(evenlist) + sequenceddata[i][tester2 + 1:]
                else:
                    if sequenceddata[i][tester2] not in oddlist:
                        sequenceddata[i] = sequenceddata[i][:tester2] + random.choice(oddlist) + sequenceddata[i][tester2 + 1:]

            indexcontainingdna = sequenceddata[i][:self.params.dna_barcode_length]
            indexbinary = dna_to_binary_custom(indexcontainingdna)
            index2list = []

            for j in range(len(indexbinary)):
                index2list.append(int(indexbinary[j]))
            index2list.reverse()
            sortablelist = []
            sortablelist.append(sequenceddata[i][self.params.dna_barcode_length:])
            for j in range(len(index2list)):
                sortablelist.insert(0, index2list[j])                    
            
            lengthofthefirst = len(index2list)
            sortedlist.append(sortablelist)
        
        listssorted = sort_lists_by_first_n_entries_synth(sortedlist, lengthofthefirst) 
        for i in range(len(listssorted)):
            for j in range(len(listssorted[i])):
                for k in range(self.params.dna_barcode_length):
                    listssorted[i][j].pop(0)
        
        listoflikley = []
        
        for evrycodeword in listssorted:
            listoflikley.append(most_common_string(evrycodeword))

        info = {}

        return listoflikley, info
