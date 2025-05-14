from dnabyte.data import DecodedData, RawData
from dnabyte.ltcodefixedsize import encode_lt, decode_lt
from reedsolo import RSCodec
from collections import defaultdict, Counter
from typing import List, Dict, Tuple
import random
import statistics
import numpy as np
import scipy as sp
from pyxdameraulevenshtein import damerau_levenshtein_distance

def dna_to_base3_no_homopolymers(dna_sequence):
    # Define the reverse mapping based on the preceding DNA base
    dna_mapping = {
        'A': {'C': '0', 'G': '1', 'T': '2'},
        'C': {'A': '0', 'G': '1', 'T': '2'},
        'G': {'A': '0', 'C': '1', 'T': '2'},
        'T': {'A': '0', 'C': '1', 'G': '2'},
        None: {'A': '0', 'C': '1', 'G': '2'}  # Initial mapping
    }
    
    base3_string = []
    previous_base = None
    
    for base in dna_sequence:
        digit = dna_mapping[previous_base][base]
        base3_string.append(digit)
        previous_base = base
    
    return ''.join(base3_string)

def base3_to_binary( base3_string):
    # Convert base-3 string to integer
    decimal_number = int(base3_string, 3)
    
    # Convert integer to binary string
    binary_string = bin(decimal_number)[2:]
    

    return binary_string

def dna_to_binary_no_homopolymers(dna_sequence):
    base3_string = dna_to_base3_no_homopolymers(dna_sequence)
    binary_string = base3_to_binary(base3_string)
    return binary_string


def binary_to_base3(binary_string):
    # Convert binary string to integer
    decimal_number = int(binary_string, 2)
    
    # Convert integer to base-3 string
    base3_string = ''
    while decimal_number > 0:
        base3_string = str(decimal_number % 3) + base3_string
        decimal_number //= 3
    
    # Ensure the base-3 string has the same length as the binary string divided by 2
    while len(base3_string) < len(binary_string) // 2:
        base3_string = '0' + base3_string
    
    return base3_string

def base3_to_dna_no_homopolymers(base3_string):
    # Define the mapping based on the preceding DNA base
    dna_mapping = {
        'A': {'0': 'C', '1': 'G', '2': 'T'},
        'C': {'0': 'A', '1': 'G', '2': 'T'},
        'G': {'0': 'A', '1': 'C', '2': 'T'},
        'T': {'0': 'A', '1': 'C', '2': 'G'},
        None: {'0': 'A', '1': 'C', '2': 'G'}  # Initial mapping
    }
    
    dna_sequence = []
    previous_base = None
    
    for digit in base3_string:
        current_base = dna_mapping[previous_base][digit]
        dna_sequence.append(current_base)
        previous_base = current_base
    
    return ''.join(dna_sequence)

def binary_to_dna_no_homopolymers(binary_string):
    base3_string = binary_to_base3(binary_string)
    dna_sequence = base3_to_dna_no_homopolymers(base3_string)
    return dna_sequence

def decimal_to_binary(n): 
    return bin(n).replace("0b", "") 

def transpose_lists(list_of_lists):
    # Use the zip function with argument unpacking to transpose the list of lists
    transposed = [list(group) for group in zip(*list_of_lists)]
    return transposed


def create_counter_list(n, m, base10_input):
    counter_list = [0] * n
    index = 0

    while base10_input > 0 and index < n:
        counter_list[index] = base10_input % m
        base10_input //= m
        index += 1

    counter_list = list(reversed(counter_list))

    return counter_list


def flatten_at_layer(nested_list, layer):
    if layer == 0:
        return nested_list
    
    flattened = []
    for element in nested_list:
        if isinstance(element, list) and layer > 0:
            flattened.extend(flatten_at_layer(element, layer - 1))
        else:
            flattened.append(element)
    
    return flattened

def bitstring_to_bytes_nobytearray(s):
    return int(s, 2).to_bytes((len(s) + 7) // 8, byteorder='big')

def bytes_to_bitstring(byte_array):
    return ''.join(format(byte, '08b') for byte in byte_array)

def bitstring_to_bytes(s):
    v = int(s, 2)
    b = bytearray()
    while v:
        b.append(v & 0xff)
        v >>= 8
    return bytes(b[::-1])


def compare_binary_strings(bin_str1, bin_str2):
    """
    Helper function for compare(). Compare two binary strings and return the indices of the differing bits.
    """
    differences = []
    min_length = min(len(bin_str1), len(bin_str2))
    
    for i in range(min_length):
        if bin_str1[i] != bin_str2[i]:
            differences.append(i)
    
    return differences


def compare(data_dec, data_raw, logger=None):
    """
    Compare decoded data with raw data and return the result of the comparison.

    Args:
        data_dec (DecodedData): The decoded data object.
        data_raw (RawData): The raw data object.

    Returns:
        str: A message indicating whether the data was successfully decoded or not.
    """
    res = compare_binary_strings(data_dec.data, data_raw.data)

    if res == [] and len(data_dec.data) == len(data_raw.data):
        return 'SUCCESS', res
    else:
        return 'ERROR', res

def compare_strings(string1, string2):
    """
    Compare two strings and return the indices where they differ.

    :param string1: The first string to compare.
    :param string2: The second string to compare.
    :return: A list of indices where the two strings differ.
    """
    # Ensure both strings are of the same length by truncating the longer one
    min_length = min(len(string1), len(string2))
    string1 = string1[:min_length]
    string2 = string2[:min_length]

    # Find the indices where the characters differ
    differences = [i for i in range(min_length) if string1[i] != string2[i]]

    return differences

def complementmap(string):
    compliment = ''
    for i in string:
        if i == 'A':
            compliment += 'T'
        if i == 'T':
            compliment += 'A'
        if i == 'C':
            compliment += 'G'
        if i == 'G':
            compliment += 'C'
    compliment = compliment[::-1]
    return compliment

def transpose_matrix(matrix):
    if not matrix:
        ValueError("Matrix is empty")
    return list(map(list, zip(*matrix)))

def split_string_into_chunks(string, n):
    # Split the string into chunks of length n
    return [string[i:i+n] for i in range(0, len(string), n)]

def makeltcode(data, num_symbols, indexcarrylength, ninputlength,howmanybitsin):
    
    makingltcode = []
    for codewords in data:
        makingltcode.append(''.join(codewords))
    
    bitsofindexcarry = indexcarrylength*howmanybitsin
    encodedltcode = encode_lt(makingltcode, len(makingltcode)+int(len(makingltcode)*num_symbols),bitsofindexcarry,ninputlength*howmanybitsin )
    splitltcode = []
    for stringies in encodedltcode:
        splitltcode.append([stringies[i:i+howmanybitsin] for i in range(0, len(stringies), howmanybitsin)])
    return splitltcode

def makeltcodesynth(data, num_symbols, indexcarrylength, ninputlength, howmanybitsin):
    bitsofindexcarry = indexcarrylength * howmanybitsin
    encodedltcode = encode_lt(data, len(data) + int(len(data) * num_symbols), bitsofindexcarry, ninputlength * howmanybitsin)
    
    return encodedltcode

def undoltcode(data, indexcarrylength, ninputlength,howmanybitsin):
    encoded_symbols = []
    for codewords in data:
        binary_list = [bin(x)[2:].zfill(howmanybitsin) for x in codewords]
        binarystring = ''.join(binary_list)
        encoded_symbols.append(binarystring)
    bitsofindexcarry = indexcarrylength*howmanybitsin
    decodedltbin, checkforvalidlt=decode_lt(encoded_symbols,bitsofindexcarry,ninputlength*howmanybitsin)
    finischeddeclist = []
     
    for decodedlt in decodedltbin:
        decodedltsplit = [decodedlt[i:i+howmanybitsin] for i in range(0, len(decodedlt), howmanybitsin)]
        decimal_list = [int(b, 2) for b in decodedltsplit]
        finischeddeclist.append(decimal_list)
    return finischeddeclist, checkforvalidlt


def undoltcodesynth(data, indexcarrylength, ninputlength,howmanybitsin):
    
    bitsofindexcarry = indexcarrylength * howmanybitsin
    decodedltbin, checkforvalidlt = decode_lt(data, bitsofindexcarry, ninputlength * howmanybitsin)
    
    
    return decodedltbin, checkforvalidlt

def bitstring_to_bytearray(bitstring):
    """
    Converts a bit string (e.g., '11001010') into a bytearray.

    Args:
        bitstring (str): A string of '0's and '1's representing bits.

    Returns:
        bytearray: The corresponding bytearray representation.
    """
    if len(bitstring) % 8 != 0:
        raise ValueError("Bit string length must be a multiple of 8")

    # Group the bit string into chunks of 8 bits and convert to bytes
    byte_array = bytearray(int(bitstring[i:i+8], 2) for i in range(0, len(bitstring), 8))
    return byte_array

def MakeReedSolomonCode(data,codewordlength,bitsinpos):
    
    reedsolomonlength = (codewordlength)*bitsinpos//8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    for codewords in data:
        bytearr = bitstring_to_bytearray(''.join(codewords))
        codeword = rs.encode(bytearr)
        reedsolobitstring = bytearray_to_bitstring(codeword)
        sploited = [reedsolobitstring[i:i+bitsinpos] for i in range(0, len(reedsolobitstring), bitsinpos)]
        reedsolomonencodedwords.append(sploited)
    
    return reedsolomonencodedwords

def MakeReedSolomonCodeSynthesis(data, errorlength):
    reedsolomonlength = (errorlength) // 8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    for codewords in data:
        bytearr = bitstring_to_bytearray(codewords)

        codewordrs = rs.encode(bytearr)
        reedsolobitstring = bytearray_to_bitstring(codewordrs)
        reedsolomonencodedwords.append(reedsolobitstring)
        
    return reedsolomonencodedwords

def bytearray_to_bitstring(byte_array):
    """
    Converts a bytearray into a bit string (e.g., '11001010').

    Args:
        byte_array (bytearray): The bytearray to convert.

    Returns:
        str: A string of '0's and '1's representing the bits.
    """
    return ''.join(format(byte, '08b') for byte in byte_array)

def undoreedsolomonsynthesis(data,errorlength):
    reedsolomonlength = (errorlength)//8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    valuechack = True
    for codewords in data:
        try:
            bytearr = bitstring_to_bytearray(codewords)
            decoded_bytearr = rs.decode(bytearr)
            decoded_bitstring = bytearray_to_bitstring(decoded_bytearr[0])
            reedsolomonencodedwords.append(decoded_bitstring)
        except:
            #logger.debug(f'Error in Reed Solomon decoding: codeword={codewords}, errorlength={errorlength}')
            valuechack = False
            reedsolomonencodedwords.append(codewords[len(codewords)-errorlength:])
    return reedsolomonencodedwords, valuechack

def undoreedsolomon(data,codewordlength,bitsinpos):
    
    reedsolomonlength = (codewordlength)*bitsinpos//8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    valuechack = True
    for codewords in data:
        # try:
        binary_list = [bin(x)[2:].zfill(bitsinpos) for x in codewords]
        bytearr = bitstring_to_bytearray(''.join(binary_list))
        decoded_bytearr = rs.decode(bytearr)
        decoded_bitstring = bytearray_to_bitstring(decoded_bytearr[0])
        sploited = [decoded_bitstring[i:i+bitsinpos] for i in range(0, len(decoded_bitstring), bitsinpos)]
        decimal_list = [int(b, 2) for b in sploited]
        reedsolomonencodedwords.append(decimal_list)
        # except:
        #     #logger.debug(f'Error in Reed Solomon decoding: codeword={codewords}, errorlength={bitsinpos}')
        #     valuechack = False
        #     reedsolomonencodedwords.append(codewords[len(codewords)-bitsinpos:])
    
    return reedsolomonencodedwords, valuechack

def check_first_n_entries_are_zero(lst, n):
    if n <= 0:
        raise ValueError("n must be greater than 0")
    
    return all(x == 0 for x in lst[:n])

def remove_complements(dna_list):
        

    unique_sequences = set()
    for seq in dna_list:
        comp_seq = complementmap(seq)
        if seq not in unique_sequences and comp_seq not in unique_sequences:
            unique_sequences.add(seq)
    
    return list(unique_sequences)

def split_string(s, n):
    first_part = s[:n]
    second_part = s[n:]
    return first_part, second_part

def binarize_and_combine_first_n_elements(lst, n, m):
    if n <= 0:
        raise ValueError("n must be greater than 0")
    
    # Extract the first n elements
    first_n_elements = lst[:n]
    
    # Convert each element to a binary string, remove the '0b' prefix, and zero-pad to length m
    binary_strings = [bin(x)[2:].zfill(m) for x in first_n_elements]
    
    # Combine the binary strings
    combined_binary_string = ''.join(binary_strings)
    
    # Convert the combined binary string back to a decimal number
    decimal_number = int(combined_binary_string, 2)
    
    return decimal_number


def count_unique_lists(list_of_lists: List[List]) -> Dict[Tuple, int]:
    # Use a defaultdict to count occurrences
    count_dict = defaultdict(int)
    
    for lst in list_of_lists:
        # Convert the list to a tuple to use it as a dictionary key
        lst_tuple = tuple(lst)
        count_dict[lst_tuple] += 1
    
    # Convert defaultdict to a regular dictionary
    return dict(count_dict)

def split_string_alternating_lengths(string, n, m):
    result = []
    i = 0
    toggle = False  # Start with n length
    while i < len(string):
        if toggle:
            substring = string[i:i+n]
            if len(substring) == n:
                result.append(substring)
            i += n
        else:
            i += m
        toggle = not toggle  # Alternate between n and m lengths
    return result


def split_string_into_chunks_poly_chain(string, ngeneric, npossition, nmessage):
    result = []
    i = 0
    length = len(string)

    # Discard the first n + m characters
    i += ngeneric + npossition
    if i >= length:
        return result

    while i < length:
        # Save the next k characters
        substring = string[i:i+nmessage]
        if len(substring) == nmessage:
            result.append(substring)
        i += nmessage
        if i >= length:
            break

        # Discard the next 2n + m characters
        i += 2 * ngeneric + npossition

    return result


def find_closest_string(target, string_list):
    min_distance = float('inf')
    closest_strings = []

    for s in string_list:
        distance = damerau_levenshtein_distance(target, s)
        if distance < min_distance:
            min_distance = distance
            closest_strings = [s]
        elif distance == min_distance:
            closest_strings.append(s)

    return random.choice(closest_strings)

def count_each_list_occurrences(list_of_lists: List[List]) -> Dict[Tuple, int]:
    """
    Counts how often each list appears in a list of lists.
    """
    count_dict = defaultdict(int)
    
    for lst in list_of_lists:
        # Convert the list to a tuple to use it as a dictionary key
        lst_tuple = tuple(lst)
        count_dict[lst_tuple] += 1
    
    # Convert defaultdict to a regular dictionary
    return dict(count_dict)

def sort_lists_by_first_n_entries_synth(lists, n):
    # Sort the entire list of lists based on the first n entries
    lists.sort(key = lambda x: x[:n])
    grouped_lists = {}
    allheaders = []
    for lst in lists:
        allheaders.append(lst[:n])
    counts = count_each_list_occurrences(allheaders)
    keys, values = zip(*counts.items())
    values = list(values)
    median_of_magnitude = sp.ndimage.median(values)
    mad_of_headers = sp.stats.median_abs_deviation(values)

    for lst in lists:
        key = tuple(lst[:n])
        if key not in grouped_lists:
            grouped_lists[key] = []
        grouped_lists[key].append(lst)

    lister = list(grouped_lists.values())
   
    listend = []
    for i in range(len(lister)):
        if len(lister[i]) >= median_of_magnitude - 5 * mad_of_headers and len(lister[i]) <= median_of_magnitude + 5 * mad_of_headers:
            listend.append(lister[i])

    return listend


def sort_lists_by_first_n_entries(lists, n, theory):
    # Sort the entire list of lists based on the first n entries
    for lst in lists:
        for k in range(n):
            if isinstance(lst[k], str):
                lists.remove(lst)
    lists.sort(key=lambda x: x[:n])
    
    grouped_lists = {}
    
    for lst in lists:
        key = tuple(lst[:n])
        if key not in grouped_lists:
            grouped_lists[key] = []
        grouped_lists[key].append(lst)

    lister = list(grouped_lists.values())
    
    if theory == 'yes':
        return lister
    else:
        listend = []
        for i in range(len(lister)):
            if len(lister[i])>0:
                listend.append(lister[i])

        return listend
    

def reduce_to_n_most_used_elements(list_of_lists, n):
    # Initialize a list to store the counts of elements at each position
    position_counts = defaultdict(Counter)
    
    # Count the occurrences of elements at each position
    for sublist in list_of_lists:
        for position, element in enumerate(sublist):
            position_counts[position][element] += 1
    
    # Initialize a list with n empty lists
    reduced_lists = [[] for _ in range(n)]
    
    # Loop through each position in position_counts
    for position in sorted(position_counts.keys()):
        # Get the n most common elements at this position
        most_common_elements = position_counts[position].most_common(n)
        if len(most_common_elements) == 1:
            for i in range(n-1):
                most_common_elements.append((most_common_elements[0][0], 1))
        
        # Distribute these elements across the n lists
        for i, (element, _) in enumerate(most_common_elements):
            reduced_lists[i].append(element)
    
    return reduced_lists


def pick_random_element_excluding(lst, exclude):
    # Filter the list to exclude the specified element
    filtered_list = [element for element in lst if element != exclude]
    
    # If the filtered list is empty, return None or raise an exception
    if not filtered_list:
        raise ValueError("No elements to choose from after excluding the specified element.")
    
    # Pick a random element from the filtered list
    return random.choice(filtered_list)


def dna_to_binary(dna_string):
    
    binary_mapping = {
        'A': '00',
        'G': '01',
        'C': '10',
        'T': '11'
    }
    binary_string = ''.join(binary_mapping[base] for base in dna_string)
    return binary_string


def bitstring_to_base4_list(bitstring):
    base4_list = [int(bitstring[i:i+2], 2) for i in range(0, len(bitstring), 2)]
    return base4_list


def flatten_single_element_lists(nested_list):
    return [inner_list[0] for inner_list in nested_list if len(inner_list) == 1]


def bitstring_to_base3_list(bitstring):
    base3_list = []
    for i in range(0, len(bitstring), 2):
        base3_digit = int(bitstring[i:i+2], 2)
        base3_list.append(base3_digit)
    return base3_list


def most_common_length(strings):
    lengths = [len(s) for s in strings]
    length_counts = Counter(lengths)
    most_common_length = length_counts.most_common(1)[0][0]
    return most_common_length



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