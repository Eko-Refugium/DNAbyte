from dnabyte.params import Params
from scipy.constants import Avogadro
import random
from tests.testlibraries.libcreation_linear_assembly import simplecration
from tests.testlibraries.libcreation_positional_assembly import generate_positional_library
import random
from simulations.auxiliary import create_one_text_file

def generate_random_integers_unsorted(total, n):
    """
    Generates n random integers such that their sum equals total, 
    each integer is at least 1, and they are randomly ordered.

    Args:
        total (int): The total sum of the integers.
        n (int): The number of integers to generate.

    Returns:
        list: A list of n integers whose sum equals total, in random order.
    """
    if n <= 0 or total < n:
        raise ValueError("n must be greater than 0 and total must be at least n.")
    
    # Adjust total to account for the minimum value of 1 for each entry
    adjusted_total = total - n
    
    # Generate n-1 random breakpoints between 0 and adjusted_total
    breakpoints = sorted(random.randint(0, adjusted_total) for _ in range(n - 1))
    
    # Add 0 and adjusted_total as the boundaries
    breakpoints = [0] + breakpoints + [adjusted_total]
    
    # Calculate the differences between consecutive breakpoints
    result = [breakpoints[i + 1] - breakpoints[i] + 1 for i in range(n)]
    
    # Shuffle the result to randomize the order
    random.shuffle(result)
    
    return result


def create_random_valid_paramset():

    assembly_structures = ['linear_assembly', 'positional_assembly']
    encoding_schemes = ['binomial_encoding', 'linear_encoding']
    outer_error_corrections = [None, 'reedsolomon']
    inner_error_corrections = [None, 'ltcode']
    sequencing_methods = [41, 40, 37, 36, 39, 38, 35,'iid', None]
    synthesis_methods = [3, 4, 5, 6, 7, 68, 69, 70, 71, 'nosynthpoly', None]
    storage_conditions = ['biogene', 'random', 'permafrost', None]

    assembly_structure_picked = random.choice(assembly_structures)
    encoding_scheme_picked = random.choice(encoding_schemes)

    if assembly_structure_picked == 'linear_assembly':
        oligo_length = random.randint(20, 100)
        motive_amount = random.randint(2, 30)
        amount_of_oligos = (2 * motive_amount) ** 2
        simplecration(oligo_length, motive_amount)
        assembly_library_name = f'lib_simple_{oligo_length}bp_{amount_of_oligos}.csv'
        if encoding_scheme_picked == 'linear_encoding':
            codeword_length_picked = random.randint(20, 160)
        elif encoding_scheme_picked == 'binomial_encoding':
            codeword_length_picked = 2 * motive_amount - 1
        else:
            raise ValueError("Invalid encoding scheme.")

    elif assembly_structure_picked == 'positional_assembly':
        length_message = random.randint(20, 50)
        amount_of_messages = random.randint(20, 250)
        length_of_generic = random.choice([i for i in range(20, 50) if i != length_message])
        length_of_position = random.choice([i for i in range(20, 50) if i != amount_of_messages and i != length_of_generic])
        amount_of_positions = random.randint(2, 100)
        generate_positional_library(length_message, amount_of_messages, length_of_generic, length_of_position, amount_of_positions)
        assembly_library_name = f'lib_positional_{length_message}({amount_of_messages})m_{length_of_generic}(2)g_{length_of_position}({amount_of_positions})p.csv'
        codeword_length_picked = amount_of_positions- 1
    else:
        raise ValueError("Invalid assembly structure.")
    
    
    outer_error_correction_picked = random.choice(outer_error_corrections)
    inner_error_correction_picked = random.choice(inner_error_corrections)

    if outer_error_correction_picked == 'reedsolomon':
        reed_solo_percentage_picked = random.uniform(0.5, 0.95)
    else:
        reed_solo_percentage_picked = None
    if inner_error_correction_picked == 'ltcode':
        amount_of_parts_of_codeword = 5
    else:
        amount_of_parts_of_codeword = 2
    
    message_carring_length = random.randint(codeword_length_picked//3, codeword_length_picked-amount_of_parts_of_codeword)
    remaningparts = codeword_length_picked - message_carring_length
    part_splits = generate_random_integers_unsorted(remaningparts, amount_of_parts_of_codeword)

    if inner_error_correction_picked != 'ltcode':
        for i in range(3):
            part_splits.append(None)

    sequencing_methods_picked = random.choice(sequencing_methods)
    if sequencing_methods_picked == 'iid':
        iid_error_rate_picked = random.uniform(0, 0.5)
    else:
        iid_error_rate_picked = None

    synthesis_methods_picked = random.choice(synthesis_methods)

    storage_conditions_picked = random.choice(storage_conditions)

    year_picked = random.randint(0, 100000)

    if encoding_scheme_picked == 'binomial_encoding':
        sigma_amount_picked = random.randint(1, 5)
    else:
        sigma_amount_picked = None

    mean_picked = random.randint(1, 100)
    std_dev_picked = random.randint(0, mean_picked//8)

    filename = create_one_text_file("./tests/testfiles/", random.randint(1, 1000))

    params_list = [Params(
        name='sim_random_params',
        filename=filename,
        assembly_structure=assembly_structure_picked, #####
        encoding_scheme=encoding_scheme_picked, #####
        library_name=assembly_library_name,  #####
        mean=mean_picked, #####
        vol=1000000 / Avogadro,
        std_dev=std_dev_picked, #####
        hybridisation_steps=10000,
        inner_error_correction=inner_error_correction_picked, #####
        outer_error_correction=outer_error_correction_picked, #####
        dna_barcode_length=part_splits[0],  #####
        codeword_maxlength_positions=part_splits[1], #####
        years=year_picked, #####
        storage_conditions=storage_conditions_picked,
        codeword_length=codeword_length_picked, #####
        percent_of_symbols=part_splits[2], #####
        ltcode_header=part_splits[3], #####
        index_carry_length=part_splits[4], #####
        synthesis_method=synthesis_methods_picked, #####
        sequencing_method=sequencing_methods_picked, #####
        iid_error_rate=iid_error_rate_picked, #####
        reed_solo_percentage=reed_solo_percentage_picked, #####
        sigma_amount=sigma_amount_picked, #####
        theory='no'
    )]

    return params_list