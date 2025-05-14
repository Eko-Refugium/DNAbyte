from pyxdameraulevenshtein import damerau_levenshtein_distance
import random
import math as m
from collections import defaultdict, Counter
from typing import List, Dict, Tuple

from dnabyte.oligo import translate_nested_list, create_positional_libraries
from dnabyte.auxiliary import complementmap, transpose_lists, transpose_matrix, split_string_into_chunks, dna_to_base3_no_homopolymers, binary_to_dna_no_homopolymers
from dnabyte.encode import Encode
from dnabyte.data import CorrectedData

class Process:

    """
    This class processes the sequenced data to data that can be decoded.
    """
    def __init__(self, library, encoding):
        self.encoding = encoding
        self.library = library


        # TODO: Check for erroneous input


    # def split_string_into_chunks(self, string, n):
    #     # Split the string into chunks of length n
    #     return [string[i:i+n] for i in range(0, len(string), n)]


    # def count_unique_lists(list_of_lists: List[List]) -> Dict[Tuple, int]:
    #     # Use a defaultdict to count occurrences
    #     count_dict = defaultdict(int)
        
    #     for lst in list_of_lists:
    #         # Convert the list to a tuple to use it as a dictionary key
    #         lst_tuple = tuple(lst)
    #         count_dict[lst_tuple] += 1
        
    #     # Convert defaultdict to a regular dictionary
    #     return dict(count_dict)


    # def split_string_alternating_lengths(self, string, n, m):
    #     result = []
    #     i = 0
    #     toggle = False  # Start with n length
    #     while i < len(string):
    #         if toggle:
    #             substring = string[i:i+n]
    #             if len(substring) == n:
    #                 result.append(substring)
    #             i += n
    #         else:
    #             i += m
    #         toggle = not toggle  # Alternate between n and m lengths
    #     return result


    # def split_string_into_chunks_poly_chain(self, string, ngeneric, npossition, nmessage):
    #     result = []
    #     i = 0
    #     length = len(string)

    #     # Discard the first n + m characters
    #     i += ngeneric + npossition
    #     if i >= length:
    #         return result

    #     while i < length:
    #         # Save the next k characters
    #         substring = string[i:i+nmessage]
    #         if len(substring) == nmessage:
    #             result.append(substring)
    #         i += nmessage
    #         if i >= length:
    #             break

    #         # Discard the next 2n + m characters
    #         i += 2 * ngeneric + npossition

    #     return result


    # def find_closest_string(self, target, string_list):
    #     min_distance = float('inf')
    #     closest_strings = []

    #     for s in string_list:
    #         distance = damerau_levenshtein_distance(target, s)
    #         if distance < min_distance:
    #             min_distance = distance
    #             closest_strings = [s]
    #         elif distance == min_distance:
    #             closest_strings.append(s)

    #     return random.choice(closest_strings)
    

    # def sort_lists_by_first_n_entries_synth(self, lists, n):
    #     # Sort the entire list of lists based on the first n entries
    #     lists.sort(key=lambda x: x[:n])
    #     grouped_lists = {}
        
    #     for lst in lists:
    #         key = tuple(lst[:n])
    #         if key not in grouped_lists:
    #             grouped_lists[key] = []
    #         grouped_lists[key].append(lst)

    #     lister = list(grouped_lists.values())
    #     listend = []
    #     for i in range(len(lister)):
    #         if len(lister[i])>0:
    #             listend.append(lister[i])

    #     return listend


    # def sort_lists_by_first_n_entries(self, lists, n, theory):
    #     # Sort the entire list of lists based on the first n entries
    #     for lst in lists:
    #         for k in range(n):
    #             if isinstance(lst[k], str):
    #                 lists.remove(lst)
    #     lists.sort(key=lambda x: x[:n])
        
    #     grouped_lists = {}
        
    #     for lst in lists:
    #         key = tuple(lst[:n])
    #         if key not in grouped_lists:
    #             grouped_lists[key] = []
    #         grouped_lists[key].append(lst)

    #     lister = list(grouped_lists.values())
        
    #     if theory == 'yes':
    #         return lister
    #     else:
    #         listend = []
    #         for i in range(len(lister)):
    #             if len(lister[i])>0:
    #                 listend.append(lister[i])

    #         return listend


    # def reduce_to_n_most_used_elements(self, list_of_lists, n):
    #     # Initialize a list to store the counts of elements at each position
    #     position_counts = defaultdict(Counter)
        
    #     # Count the occurrences of elements at each position
    #     for sublist in list_of_lists:
    #         for position, element in enumerate(sublist):
    #             position_counts[position][element] += 1
        
    #     # Initialize a list with n empty lists
    #     reduced_lists = [[] for _ in range(n)]
        
    #     # Loop through each position in position_counts
    #     for position in sorted(position_counts.keys()):
    #         # Get the n most common elements at this position
    #         most_common_elements = position_counts[position].most_common(n)
            
    #         # Distribute these elements across the n lists
    #         for i, (element, _) in enumerate(most_common_elements):
    #             reduced_lists[i].append(element)
        
    #     return reduced_lists


    # def pick_random_element_excluding(self, lst, exclude):
    #     # Filter the list to exclude the specified element
    #     filtered_list = [element for element in lst if element != exclude]
        
    #     # If the filtered list is empty, return None or raise an exception
    #     if not filtered_list:
    #         raise ValueError("No elements to choose from after excluding the specified element.")
        
    #     # Pick a random element from the filtered list
    #     return random.choice(filtered_list)


    # def dna_to_binary(self, dna_string):
        
    #     binary_mapping = {
    #         'A': '00',
    #         'G': '01',
    #         'C': '10',
    #         'T': '11'
    #     }
    #     binary_string = ''.join(binary_mapping[base] for base in dna_string)
    #     return binary_string
    
    # def bitstring_to_base4_list(self, bitstring):
    #     base4_list = [int(bitstring[i:i+2], 2) for i in range(0, len(bitstring), 2)]
    #     return base4_list
    

    # def flatten_single_element_lists(self, nested_list):
    #     return [inner_list[0] for inner_list in nested_list if len(inner_list) == 1]
    

    # def bitstring_to_base3_list(self, bitstring):
    #     base3_list = []
    #     for i in range(0, len(bitstring), 2):
    #         base3_digit = int(bitstring[i:i+2], 2)
    #         base3_list.append(base3_digit)
    #     return base3_list


    # def most_common_length(self, strings):
    #     lengths = [len(s) for s in strings]
    #     length_counts = Counter(lengths)
    #     most_common_length = length_counts.most_common(1)[0][0]
    #     return most_common_length
    

    # def most_common_string(self, strings):
    #     if not strings:
    #         return ""
        
    #     # Initialize a list to store the most common letters at each position
    #     most_common_letters = []
        
    #     # Iterate over each position in the strings
    #     for i in range(len(strings[0])):
    #         # Get the letters at the current position from all strings
    #         letters_at_position = [s[i] for s in strings]
            
    #         # Find the most common letter at the current position
    #         most_common_letter = Counter(letters_at_position).most_common(1)[0][0]
            
    #         # Append the most common letter to the list
    #         most_common_letters.append(most_common_letter)
        
    #     # Join the list of most common letters into a single string
    #     return ''.join(most_common_letters)


    # def dna_to_binary_custom(self, dna_string):
    #     binary_string = []
    #     for i, base in enumerate(dna_string):
    #         if i % 2 == 0:  # Even position
    #             if base == 'T':
    #                 binary_string.append('0')
    #             elif base == 'G':
    #                 binary_string.append('1')
    #         else:  # Odd position
    #             if base == 'C':
    #                 binary_string.append('0')
    #             elif base == 'A':
    #                 binary_string.append('1')
    #     return ''.join(binary_string)


    def process(self, sequenceddata, barcodelength, **kwargs):

        sequenceddata = sequenceddata.data
        library = self.library
        encodingscheme = self.encoding


        if library == 'synthesis':
            if encodingscheme == 'max_density_encoding':

                if kwargs['outer_error_correction'] == 'reedsolomon':
                    barcode = barcodelength
                    # while barcode%4 != 0:
                    #     barcode += 1
                    tester = kwargs['codewordlength']
                    # while tester%4 != 0:
                    #     tester += 1
                
                else:
                    barcode = barcodelength
                    tester = kwargs['codewordlength']
                sortedlist = []
                sequenceddata = [seq.replace(' ', '') for seq in sequenceddata] 
                

                # TODO: parallelize
                for i in range(len(sequenceddata)):
                    
                    bases = ['A', 'C', 'T', 'G']
                    while len(sequenceddata[i]) < tester:
                        # Randomly choose a position to insert a base
                        insert_position = random.randint(0, len(sequenceddata[i]))
                        # Randomly choose a base to insert
                        base_to_insert = random.choice(bases)
                        # Insert the base at the chosen position
                        sequenceddata[i] = sequenceddata[i][:insert_position] + base_to_insert + sequenceddata[i][insert_position:]
                    
                    while len(sequenceddata[i]) > tester:
                        # Randomly choose a position to delete a base
                        delete_position = random.randint(0, len(sequenceddata[i]) - 1)
                        # Delete the base at the chosen position
                        sequenceddata[i] = sequenceddata[i][:delete_position] + sequenceddata[i][delete_position + 1:]
                    
                    indexcontainingdna = sequenceddata[i][:barcode]
                    indexbinary = self.dna_to_binary(indexcontainingdna)
                    indexb4list = self.bitstring_to_base4_list(indexbinary)
                    indexb4list.reverse()
                    sortablelist = []
                    sortablelist.append(sequenceddata[i][barcode:])
                    for j in range(len(indexb4list)):
                        sortablelist.insert(0, indexb4list[j])

                    lengthofthefirst = len(indexb4list)
                    sortedlist.append(sortablelist)

                listssorted = self.sort_lists_by_first_n_entries_synth(sortedlist, lengthofthefirst) 

                for i in range(len(listssorted)):
                    for j in range(len(listssorted[i])):
                        for k in range(barcode):
                            listssorted[i][j].pop(0)
                
                listoflikley = []
                
                
                for evrycodeword in listssorted:
                    listoflikley.append(self.most_common_string(evrycodeword))
                # for i in range(len(listoflikley)):
                #     listoflikley[i] = listoflikley[i][0]
                # # listoflikley = self.flatten_single_element_lists(listoflikley)
                
                info = {}

                return listoflikley, info
                
            # TODO: parallelize
            elif encodingscheme == 'no_homopolymeroddeven_encoding':

                if kwargs['outer_error_correction'] == 'reedsolomon':
                    barcode = barcodelength
                    # while barcode%8 != 0:
                    #     barcode -= 1
                    tester = kwargs['codewordlength']
                    # while tester%8 != 0:
                    #     tester += 1
                
                else:
                    barcode = barcodelength
                    tester = kwargs['codewordlength']
                sortedlist = []

                sequenceddata = [seq.replace(' ', '') for seq in sequenceddata] 
                for i in range(len(sequenceddata)):
                    
                    bases = ['A', 'C', 'T', 'G']
                    
                    while len(sequenceddata[i]) < tester:
                        for tester2 in range(len(sequenceddata[i])-1):
                            evenlist = ['T', 'G']
                            oddlist = ['C', 'A']
                            if sequenceddata[i][tester2] in evenlist and sequenceddata[i][tester2+1] in evenlist:
                                sequenceddata[i] = sequenceddata[i][:tester2+1] + random.choice(oddlist) + sequenceddata[i][tester2+2:]
                                break
                            elif sequenceddata[i][tester2] in oddlist and sequenceddata[i][tester2+1] in oddlist:
                                sequenceddata[i] = sequenceddata[i][:tester2+1] + random.choice(evenlist) + sequenceddata[i][tester2+2:]
                                break
                        else:
                            # Randomly choose a position to insert a base
                            insert_position = random.randint(0, len(sequenceddata[i]))
                            # Randomly choose a base to insert
                            base_to_insert = random.choice(bases)
                            # Insert the base at the chosen position
                            sequenceddata[i] = sequenceddata[i][:insert_position] + base_to_insert + sequenceddata[i][insert_position:]
                            
                    
                    while len(sequenceddata[i]) > tester:
                        # Randomly choose a position to delete a base
                        delete_position = random.randint(0, len(sequenceddata[i]) - 1)
                        # Delete the base at the chosen position
                        sequenceddata[i] = sequenceddata[i][:delete_position] + sequenceddata[i][delete_position + 1:]
                    for tester2 in range(len(sequenceddata[i])-1):
                        evenlist = ['T', 'G']
                        oddlist = ['C', 'A']
                        if tester2%2 == 0:
                            if sequenceddata[i][tester2] not in evenlist:
                                sequenceddata[i] = sequenceddata[i][:tester2] + random.choice(evenlist) + sequenceddata[i][tester2+1:]
                        else:
                            if sequenceddata[i][tester2] not in oddlist:
                                sequenceddata[i] = sequenceddata[i][:tester2] + random.choice(oddlist) + sequenceddata[i][tester2+1:]

                    indexcontainingdna = sequenceddata[i][:barcode]
                    indexbinary = self.dna_to_binary_custom(indexcontainingdna)
                    index2list = []
                    for j in range(len(indexbinary)):
                        index2list.append(int(indexbinary[j]))
                    index2list.reverse()
                    sortablelist = []
                    sortablelist.append(sequenceddata[i][barcode:])
                    for j in range(len(index2list)):
                        sortablelist.insert(0, index2list[j])                    
                    
                    lengthofthefirst = len(index2list)
                    sortedlist.append(sortablelist)
                
                listssorted = self.sort_lists_by_first_n_entries_synth(sortedlist, lengthofthefirst) 
                for i in range(len(listssorted)):
                    for j in range(len(listssorted[i])):
                        for k in range(barcode):
                            listssorted[i][j].pop(0)
                
                listoflikley = []
                
                for evrycodeword in listssorted:
                    listoflikley.append(self.most_common_string(evrycodeword))

                info = {}

                return listoflikley, info


            elif encodingscheme == 'no_homopolymer_encoding':
                if kwargs['outer_error_correction'] == 'reedsolomon':
                    barcode = barcodelength
                    while barcode%4 != 0:
                        barcode -= 1
                    tester = kwargs['codewordlength']
                    while tester%4 != 0:
                        tester += 1
                else:
                    barcode = barcodelength
                    tester = kwargs['codewordlength']
                sortedlist = []
                sequenceddata = [seq.replace(' ', '') for seq in sequenceddata] 


                # TODO: parallize
                for i in range(len(sequenceddata)):
                    
                    bases = ['A', 'C', 'T', 'G']
                    
                    while len(sequenceddata[i]) < tester:
                        # Randomly choose a position to insert a base
                        insert_position = random.randint(0, len(sequenceddata[i]))
                        # Randomly choose a base to insert
                        base_to_insert = random.choice(bases)
                        # Insert the base at the chosen position
                        sequenceddata[i] = sequenceddata[i][:insert_position] + base_to_insert + sequenceddata[i][insert_position:]
                    
                    while len(sequenceddata[i]) > tester:
                        # Randomly choose a position to delete a base
                        delete_position = random.randint(0, len(sequenceddata[i]) - 1)
                        # Delete the base at the chosen position
                        sequenceddata[i] = sequenceddata[i][:delete_position] + sequenceddata[i][delete_position + 1:]                            
                             
                    indexcontainingdna = sequenceddata[i][:barcode]
                    for tester2 in range(len(indexcontainingdna)-1):
                        if indexcontainingdna[tester2] == indexcontainingdna[tester2+1]:
                            lst = ['A', 'C', 'T', 'G']
                            exclude = indexcontainingdna[tester2+1]
                            indexcontainingdna = indexcontainingdna[:tester2+1] + self.pick_random_element_excluding(lst, exclude) + indexcontainingdna[tester2+ 2:]

                    indexbinary = dna_to_base3_no_homopolymers(indexcontainingdna)
                    indexb3list = []
                    for j in range(len(indexbinary)):
                        indexb3list.append(int(indexbinary[j]))
                    indexb3list.reverse()
                    sortablelist = []
                    sortablelist.append(sequenceddata[i][barcode:])
                    for j in range(len(indexb3list)):
                        sortablelist.insert(0, indexb3list[j])
                    
                    lengthofthefirst = len(indexb3list)
                    sortedlist.append(sortablelist)
                
                listssorted = self.sort_lists_by_first_n_entries_synth(sortedlist, lengthofthefirst) 
                for i in range(len(listssorted)):
                    for j in range(len(listssorted[i])):
                        for k in range(barcode):
                            listssorted[i][j].pop(0)
                
                listoflikley = []
                for evrycodeword in listssorted:
                    listoflikley.append(self.most_common_string(evrycodeword))

                for i in range(len(listoflikley)-1):
                    if listoflikley[i] == listoflikley[i+1]:
                        lst = ['A', 'C', 'T', 'G']
                        exclude = listoflikley[i+1]
                        listoflikley = listoflikley[:i+1] + self.pick_random_element_excluding(lst, exclude) + listoflikley[i+ 2:]

                info = {}

                return listoflikley, info


        elif library.structure == 'linear_assembly':
            if  encodingscheme== 'linear_encoding':
                finishcedlisdtofdecodedDNA = []
                DNAs = library.library
                oligolength = len(DNAs[0])
                seperated = []
                

                # TODO: parallelize
                for i in range(len(sequenceddata)):
                    fivetothreeend = sequenceddata[i][1]
                    bases = ['A', 'C', 'T', 'G']

                    if kwargs['codewordlength']%2==0:
                        lengthofcorrect = (kwargs['codewordlength'])*(len(DNAs[0]))
                    else:
                        lengthofcorrect = kwargs['codewordlength']*(len(DNAs[0]))
                    
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
                        
                    chunks = split_string_into_chunks(fivetothreeend,oligolength)
                    correctedlistofinfooligos = []
                    for singularinformationoligo in chunks:
                        correctedlistofinfooligos.append(self.find_closest_string(singularinformationoligo, DNAs))
                    for j in range(barcodelength):
                        correctedlistofinfooligos[j]=DNAs.index(correctedlistofinfooligos[j])
                    seperated.append(correctedlistofinfooligos)
                    
                listssorted = self.sort_lists_by_first_n_entries(seperated, barcodelength,theory=kwargs['theory']) 
                for i in range(len(listssorted)):
                    for j in range(len(listssorted[i])):
                        for k in range(barcodelength):
                            listssorted[i][j].pop(0)
                listoflikley = []
                for evrycodeword in listssorted:
                    listoflikley.append(self.reduce_to_n_most_used_elements(evrycodeword,1))
                for i in range(len(listoflikley)):
                    listoflikley[i] = listoflikley[i][0]
                
                info = {}       

                return listoflikley, info
            

            if encodingscheme == 'binomial_encoding':
                finishcedlisdtofdecodedDNA = []
                DNAs = library.library
                oligolength = len(DNAs[0])
                positionallib=create_positional_libraries(library, len(library.leftmotives)+len(library.rightmotives) - 1)

                transladictleft = library.translationlibleft
                transladictright = library.translationlibright

                seperated = []

                # TODO: parallelize
                for i in range(len(sequenceddata)):
                    fivetotreeprimetuples = []
                    #for singular_codeword in codewords_tuples:
                    fivetothreeend = sequenceddata[i][0]
                    
                    lengthofcorrect = (len(library.leftmotives)+2*len(library.rightmotives)-1)*oligolength 
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
                      
                    listofinfooligos = self.split_string_alternating_lengths(fivetothreeend, oligolength, oligolength//2)
                    
                    for j in range(len(listofinfooligos)):
                        if j%2==1:
                            listofinfooligos[j] = (listofinfooligos[j][:len(listofinfooligos[j])//2][::-1]+listofinfooligos[j][len(listofinfooligos[j])//2:][::-1])[::-1]
                        else:
                            listofinfooligos[j] = complementmap(listofinfooligos[j][:len(listofinfooligos[j])//2])+complementmap(listofinfooligos[j][len(listofinfooligos[j])//2:])

                    correctedlistofinfooligos = []
                    for j, singularinformationoligo in enumerate(listofinfooligos):
                        correctedlistofinfooligos.append(self.find_closest_string(singularinformationoligo, positionallib[j]))

                    for j in range(barcodelength):
                        if correctedlistofinfooligos[j] in positionallib[j]:
                            correctedlistofinfooligos[j]=positionallib[j].index(correctedlistofinfooligos[j])
                        # else:
                        #     correctedlistofinfooligos[j]=random.randint(0,len(positionallib[j])-1)
                    seperated.append(correctedlistofinfooligos)
                listssorted = self.sort_lists_by_first_n_entries(seperated, barcodelength,theory=kwargs['theory']) 
                
                for i in range(len(listssorted)):
                    for j in range(len(listssorted[i])):
                        for k in range(barcodelength):
                            listssorted[i][j].pop(0)
                
                listoflikley = []
                for evrycodeword in listssorted:
                    listoflikley.append(self.reduce_to_n_most_used_elements(evrycodeword, kwargs['sigmaamount']))

                transposeinfo = []
                for i in range(len(listoflikley)):
                    transposeinfo.append(transpose_matrix(listoflikley[i]))
                
                info = {}

                return transposeinfo, info

            
        if library.structure == 'positional_assembly':
            if encodingscheme == 'linear_encoding':
                finishcedlisdtofdecodedDNA = []
                DNAs = library.messages
                generic = library.generic
                psitoion = library.position
                DNAs.sort()
                
                oligolength = len(DNAs[0])
                genericlength = len(generic[0])
                psitoionlength = len(psitoion[0])
                
                seperated = []
                for i in range(len(sequenceddata)):
                    
                    fivetothreeend = sequenceddata[i][0]

                    if len(psitoion)%2==0:
                        lengthofcorrect = (len(psitoion)-1)*(oligolength+2*genericlength)+(len(psitoion))*psitoionlength
                    else:
                        lengthofcorrect = (len(psitoion)-1)*(oligolength+2*genericlength)+(len(psitoion)-1)*psitoionlength
                    bases = ['A', 'C', 'T', 'G']
                    
                    while len(fivetothreeend) < lengthofcorrect:
                        # Randomly choose a position to insert a base
                        insert_position = random.randint(barcodelength, len(fivetothreeend))
                        # Randomly choose a base to insert
                        base_to_insert = random.choice(bases)
                        # Insert the base at the chosen position
                        fivetothreeend = fivetothreeend[:insert_position] + base_to_insert + fivetothreeend[insert_position:]
                    
                    while len(fivetothreeend) > lengthofcorrect:
                        # Randomly choose a position to delete a base
                        delete_position = random.randint(0, len(fivetothreeend) - 1)
                        # Delete the base at the chosen position
                        fivetothreeend = fivetothreeend[:delete_position] + fivetothreeend[delete_position + 1:]
                    
                    chunks = self.split_string_into_chunks_poly_chain(fivetothreeend,genericlength,psitoionlength,oligolength)
                    correctedlistofinfooligos = []
                    for singularinformationoligo in chunks:
                        correctedlistofinfooligos.append(self.find_closest_string(singularinformationoligo, DNAs))
                    for j in range(barcodelength):
                        correctedlistofinfooligos[j]=DNAs.index(correctedlistofinfooligos[j])
                    seperated.append(correctedlistofinfooligos)
                
                    
                listssorted = self.sort_lists_by_first_n_entries(seperated, barcodelength,theory=kwargs['theory']) 

                for i in range(len(listssorted)):
                    for j in range(len(listssorted[i])):
                        for k in range(barcodelength):
                            listssorted[i][j].pop(0)

                listoflikley = []
                for evrycodeword in listssorted:
                    listoflikley.append(self.reduce_to_n_most_used_elements(evrycodeword,1))
                for i in range(len(listoflikley)):
                    listoflikley[i] = listoflikley[i][0]
                
                info = {}

                return listoflikley, info
            
            if encodingscheme == 'binomial_encoding':
                finishcedlisdtofdecodedDNA = []
                DNAs = library.messages
                generic = library.generic
                psitoion = library.position
                oligolength = len(DNAs[0])
                genericlength = len(generic[0])
                psitoionlength = len(psitoion[0])
                DNAs.sort()
                

                seperated = []
                for i in range(len(sequenceddata)):
                    #codwordtemp = []
                    #for codewords_tuples in sequenceddata[i]:
                    fivetotreeprimetuples = []
                    #for singular_codeword in codewords_tuples:
                    fivetothreeend = sequenceddata[i][0]

                    if len(psitoion)%2==0:
                        lengthofcorrect = (len(psitoion)-1)*(oligolength+2*genericlength)+(len(psitoion))*psitoionlength
                    else:
                        lengthofcorrect = (len(psitoion)-1)*(oligolength+2*genericlength)+(len(psitoion)-1)*psitoionlength
                    bases = ['A', 'C', 'T', 'G']
                    
                    while len(fivetothreeend) < lengthofcorrect:
                        # Randomly choose a position to insert a base
                        insert_position = random.randint(barcodelength, len(fivetothreeend))
                        # Randomly choose a base to insert
                        base_to_insert = random.choice(bases)
                        # Insert the base at the chosen position
                        fivetothreeend = fivetothreeend[:insert_position] + base_to_insert + fivetothreeend[insert_position:]
                    
                    while len(fivetothreeend) > lengthofcorrect:
                        # Randomly choose a position to delete a base
                        delete_position = random.randint(0, len(fivetothreeend) - 1)
                        # Delete the base at the chosen position
                        fivetothreeend = fivetothreeend[:delete_position] + fivetothreeend[delete_position + 1:]
                    
                    
                    listofinfooligos = self.split_string_into_chunks_poly_chain(fivetothreeend,genericlength,psitoionlength,oligolength)
                    correctedlistofinfooligos = []
                    for singularinformationoligo in listofinfooligos:
                        correctedlistofinfooligos.append(self.find_closest_string(singularinformationoligo, DNAs))
                        #fivetotreeprimetuples.append(correctedlistofinfooligos)
                    for j in range(barcodelength):
                        correctedlistofinfooligos[j]=DNAs.index(correctedlistofinfooligos[j])
                    seperated.append(correctedlistofinfooligos)
                
                listssorted = self.sort_lists_by_first_n_entries(seperated, barcodelength,theory=kwargs['theory']) 
                for i in range(len(listssorted)):
                    for j in range(len(listssorted[i])):
                        for k in range(barcodelength):
                            listssorted[i][j].pop(0)
                
                
                listoflikley = []
                for evrycodeword in listssorted:

                    listoflikley.append(self.reduce_to_n_most_used_elements(evrycodeword, kwargs['sigmaamount']))

                transposeinfo = []
                for i in range(len(listoflikley)):
                    transposeinfo.append(transpose_matrix(listoflikley[i]))


                info = {}

                return transposeinfo, info
