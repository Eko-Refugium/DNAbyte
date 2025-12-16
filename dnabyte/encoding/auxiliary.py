import os
import json
from collections import defaultdict, Counter
from typing import List, Dict, Tuple
import scipy as sp
from pyxdameraulevenshtein import damerau_levenshtein_distance
import random
from typing import List, Dict, Tuple

from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.library import Library
from dnabyte.encode import Encode


def create_counter_list(n, m, base10_input):
    counter_list = [0] * n
    index = 0

    while base10_input > 0 and index < n:
        counter_list[index] = base10_input % m
        base10_input //= m
        index += 1

    counter_list = list(reversed(counter_list))

    return counter_list

def decimal_to_binary(n): 
    return bin(n).replace("0b", "") 

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

def data_as_indixes(data, library_lst):
    """
    Convert a nested list of library elements to a nested list of indices of the sequences in the library.

    Parameters:
    data (list): A nested list of library elements.
    library (list): A list of sequences.

    Returns:
    list: A nested list where the elements are the indices of the sequences in the library.
    """
    def get_index(sequence):
        try:
            return library_lst.index(sequence)
        except ValueError:
            return -1  # Return -1 if the sequence is not found in the library

    def convert_nested_list(nested_list):
        if isinstance(nested_list, list):
            return [convert_nested_list(item) for item in nested_list]
        else:
            return get_index(nested_list)

    return convert_nested_list(data)


def indices_as_data(data, library):
    """
    Convert a nested list of indices to a nested list of library elements.

    Parameters:
    data (list): A nested list of indices.
    library (list): A list of sequences.

    Returns:
    list: A nested list where the elements are the sequences from the library.
    """
    def get_sequence(index):
        try:
            return library[index]
        except IndexError:
            return None  # Return None if the index is out of range

    def convert_nested_list(nested_list):
        if isinstance(nested_list, list):
            return [convert_nested_list(item) for item in nested_list]
        else:
            return get_sequence(nested_list)

    return convert_nested_list(data)


def create_json_file(type, filenames, assembly_structure, encoding_scheme, data, 
                     inner_error_correction, outer_error_correction, 
                     theory, **kwargs):

    encoding_parameters = {
        "codeword_length": kwargs.get('codewordlength'),
        "dna_barcode_length": kwargs.get('dna_barcode_length'),
        "codeword_maxlength_positions": kwargs.get('codeword_maxlength_positions'),
        "percent_of_symbols": kwargs.get('percentofsymbols'),
        "index_carry_length": kwargs.get('indexcarrylength'),
        "theory": theory
    }

    error_correction_parameters = {
        "inner_error_correction": inner_error_correction,
        "outer_error_correction": outer_error_correction,
        "reed_solomon_percentage": kwargs.get('reedsolopercentage')
    }

    json_data = {
        "type": type,
        "filenames": filenames,
        "assembly_structure": assembly_structure,
        "encoding_scheme": encoding_scheme,
        "encoding_parameters": encoding_parameters,
        "error_correction_parameters": error_correction_parameters
    }


    if type == 'simulated_data':
        assembly_parameters = {
            "mean_molar": kwargs.get('mean_molar'),
            "std_dev_molar": kwargs.get('std_dev_molar'),
            "hybridisation_steps": kwargs.get('hybridisation_steps')
        }
        
        if assembly_structure == 'synthesis':
            assembly_parameters['synthesis_method'] = kwargs['synthesis_method']

        storage_parameters = {
            "years": kwargs.get('years'),
            "storage_conditions": kwargs.get('storage_conditions')
        }

        sequencing_parameters = {
            "sequencing_technology": kwargs.get('sequencing_technology')
        }

        json_data["assembly_parameters"] = assembly_parameters
        json_data["storage_parameters"] = storage_parameters
        json_data["sequencing_parameters"] = sequencing_parameters


    # save assembled data as indices of library elements
    if assembly_structure == 'synthesis':
        json_data['data'] = data.data
    elif assembly_structure != 'synthesis' and type != 'encoded_data':
        json_data['library_name'] = kwargs.get('library_name')
        json_data['library'] = kwargs.get('library').library
        json_data['data'] = data.data
    else:
        json_data['library_name'] = kwargs.get('library_name')
        json_data['library'] = kwargs.get('library').library
        json_data['data'] = data_as_indixes(data.data, kwargs.get('library').library)

    return json_data


# def save_json_file(data_json, directory, filename):
#     """
#     Save a JSON object to a file.

#     Parameters:
#     data_json (dict): A JSON object.
#     directory (str): The directory where the file will be saved.
#     filename (str): The name of the file.

#     Returns:
#     None
#     """
#     # Serialize parts of the JSON object with indents
#     json_part_with_indents = json.dumps({
#     "filenames": data_json["filenames"],
#     "type": data_json["type"],
#     "assembly_structure": data_json["assembly_structure"],
#     "encoding_scheme": data_json["encoding_scheme"],
#     "encoding_parameters": data_json["encoding_parameters"],
#     "error_correction_parameters": data_json["error_correction_parameters"]
#     }, indent=4)

#     # Add assembly, storage, and sequencing parameters if present
#     if data_json['type'] == 'simulated_data':
#         json_part_with_indents = json_part_with_indents.rstrip('} \n') + ',\n    ' + json.dumps({
#             "assembly_parameters": data_json["assembly_parameters"],
#             "storage_parameters": data_json["storage_parameters"],
#             "sequencing_parameters": data_json["sequencing_parameters"]
#         }, indent=4).lstrip('{')

#     # Serialize library and data of the JSON object without indents
#     if data_json['assembly_structure'] == 'synthesis':
#         json_part_without_indents = json.dumps({
#         "data": data_json["data"]
#         }, indent=None, separators=(',', ':'))
#     else:
#         json_part_without_indents = json.dumps({
#             "library": data_json["library"],
#             "data": data_json["data"]
#         }, indent=None, separators=(',', ':'))

#     # Combine the two parts
#     combined_json = json_part_with_indents.rstrip('} \n') + ',\n    ' + json_part_without_indents.lstrip('{') + '\n}'

#     # Save the combined JSON string to a file
#     json_filename = os.path.join(directory, filename)
#     with open(json_filename, 'w') as json_file:
#         json_file.write(combined_json)

# def save_json_file(data_json, directory, filename):
#     """
#     Save a JSON object to a file.

#     Parameters:
#     data_json (dict): A JSON object.
#     directory (str): The directory where the file will be saved.
#     filename (str): The name of the file.

#     Returns:
#     None
#     """
#     # Ensure the directory exists
#     if not os.path.exists(directory):
#         os.makedirs(directory)

#     # Serialize the entire JSON object with indents
#     json_string = json.dumps(data_json, indent=4)

#     # Save the JSON string to a file
#     json_filename = os.path.join(directory, filename)
#     with open(json_filename, 'w') as json_file:
#         json_file.write(json_string)


def save_json_file(data_json, directory, filename):
    """
    Save a JSON object to a file.

    Parameters:
    data_json (dict): A JSON object.
    directory (str): The directory where the file will be saved.
    filename (str): The name of the file.

    Returns:
    None
    """
    # Ensure the directory exists
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Serialize parts of the JSON object with indents
    json_part_with_indents = json.dumps({
        "filenames": data_json["filenames"],
        "type": data_json["type"],
        "assembly_structure": data_json["assembly_structure"],
        "encoding_scheme": data_json["encoding_scheme"],
        "encoding_parameters": data_json["encoding_parameters"],
        "error_correction_parameters": data_json["error_correction_parameters"],
        "assembly_parameters": data_json.get("assembly_parameters"),
        "storage_parameters": data_json.get("storage_parameters"),
        "sequencing_parameters": data_json.get("sequencing_parameters"),
        "library_name": data_json.get("library_name"),
        "library": data_json.get("library")
    }, indent=4)

    # Serialize library and data of the JSON object without indents
    json_part_without_indents = json.dumps({
        "data": data_json["data"]
    }, indent=None, separators=(',', ': '))

    # Combine the two parts
    combined_json = json_part_with_indents[:-2] + ', \n    ' + json_part_without_indents[1:] #+ '\n}'
    #combined_json = json_part_with_indents.rstrip('} \n') + '},\n    ' + json_part_without_indents.lstrip('{') #+ '\n}'
    # Save the combined JSON string to a file
    json_filename = os.path.join(directory, filename)
    with open(json_filename, 'w') as json_file:
        json_file.write(combined_json)


def create_encodeddata_object(json_data):
    """
    Create an EncodedData object from a JSON object.

    Parameters:
    json_data (dict): A JSON object.

    Returns:
    EncodedData: An EncodedData object.
    """
    if json_data['assembly_structure'] == 'synthesis':
        return NucleobaseCode(data = json_data['data'], file_paths = json_data['filenames'])

    else:
        return NucleobaseCode(data = indices_as_data(json_data['data'], json_data['library']), file_paths = json_data['filenames'])



def create_encode_object(json_data):
    """
    Create an Encode object from a JSON object.

    Parameters:
    json_data (dict): A JSON object.

    Returns:
    Encode: An Encode object.
    """

    library_name = json_data['library_name']
    # load library
    wd = os.getcwd()

    # load library
    library_path = os.path.join(wd, 'app', 'static', 'libraries', library_name + '.csv')
    lib = Library(structure=json_data['assembly_structure'], filename = library_path)
    library = Library(structure=json_data['assembly_structure'], library=lib)

    # create Encode object
    return Encode(assembly_structure = json_data['assembly_structure'], 
                  encoding_scheme = json_data['encoding_scheme'], 
                  inner_error_correction = json_data['error_correction_parameters']['inner_error_correction'],
                  outer_error_correction = json_data['error_correction_parameters']['outer_error_correction'],
                  dna_barcode_length = json_data['encoding_parameters']['dna_barcode_length'],
                  codeword_maxlength_positions = json_data['encoding_parameters']['codeword_maxlength_positions'],
                  codewordlength = json_data['encoding_parameters']['codeword_length'],
                  library = library,
                  library_name = library_name,
                  sigmaamount = json_data['encoding_parameters']['index_carry_length'],
                  theory = json_data['encoding_parameters']['theory']
                  )




def create_sequenceddata_object(json_data):

    if json_data['type'] == 'sequenced_data':
        return InSilicoDNA(data = json_data['data'])

    else:
        raise ValueError('The JSON object does not contain sequenced data.')




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

def sort_lists_by_first_n_entries(lists, n, theory):
    # Filter out lists that have strings in the first n positions (keep only numeric entries)
    filtered_lists = []
    for lst in lists:
        has_string = False
        for k in range(min(n, len(lst))):
            if isinstance(lst[k], str):
                has_string = True
                break
        if not has_string:
            filtered_lists.append(lst)
    
    # Sort the filtered lists by the first n numeric entries
    filtered_lists.sort(key=lambda x: x[:n])
    
    # Group lists that have the same first n entries
    grouped_lists = {}
    
    for lst in filtered_lists:
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


def bytes_to_bitstring(byte_array):
    return ''.join(format(byte, '08b') for byte in byte_array)






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



def transpose_matrix(matrix):
    if not matrix:
        ValueError("Matrix is empty")
    return list(map(list, zip(*matrix)))

def split_string_into_chunks(string, n):
    # Split the string into chunks of length n
    return [string[i:i+n] for i in range(0, len(string), n)]




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


def check_library(inputparams, default, assembly_structure):
    print('current dir:', os.getcwd())
    if not hasattr(inputparams, 'library_name') or inputparams.library_name is not None:
        # set default library
        library_name = default
        library = Library(structure=assembly_structure, filename='./tests/testlibraries/' + default)
    else:
        library_name = inputparams.library_name
        with open('./tests/testlibraries/' + inputparams.library_name, 'r') as f:
            first_line = f.readline().strip() 
            if first_line == 'Messages':
                raise ValueError("Library not compatible with linear assembly")
            else:
                library = Library(structure=assembly_structure, filename='./tests/testlibraries/' + inputparams.library_name)

    return library_name, library

def check_parameter(parameter, default, min, max, inputparams):
    if not hasattr(inputparams, parameter) or inputparams.__dict__[parameter] is None:
        parameter_value = default
    elif not (min <= inputparams.__dict__[parameter] <= max):
        raise ValueError(f"{parameter} must be greater than or equal to {min} and less than or equal to {max}, got {inputparams.__dict__[parameter]}")
    else:
        parameter_value = inputparams.__dict__[parameter]
    
    return parameter_value





    









