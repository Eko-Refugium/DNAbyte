import os
import json
from dnabyte.data import EncodedData, SequencedData
from dnabyte.library import Library
from dnabyte.encode import Encode


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
        "codeword_max_lenght_positions": kwargs.get('codeword_maxlength_positions'),
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
        return EncodedData(raw_data = json_data['data'], file_paths = json_data['filenames'])

    else:
        return EncodedData(raw_data = indices_as_data(json_data['data'], json_data['library']), file_paths = json_data['filenames'])



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
                  codeword_maxlength_positions = json_data['encoding_parameters']['codeword_max_lenght_positions'],
                  codewordlength = json_data['encoding_parameters']['codeword_length'],
                  library = library,
                  library_name = library_name,
                  sigmaamount = json_data['encoding_parameters']['index_carry_length'],
                  theory = json_data['encoding_parameters']['theory']
                  )




def create_sequenceddata_object(json_data):

    if json_data['type'] == 'sequenced_data':
        return SequencedData(raw_data = json_data['data'])

    else:
        raise ValueError('The JSON object does not contain sequenced data.')





