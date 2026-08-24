import os
import string
import random

from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.encode import Encode
from dnabyte.params import Params
from dnabyte.sequence import SimulateSequencing
from dnabyte.store import SimulateStorage
from dnabyte.synthesize import SimulateSynthesis
from simulations.simulation import Simulation
from dnabyte.binarize import Binarize
from dnabyte.data_classes import Data


#Function that makes a random n bit txt file
def create_random_file(output_path: str, size_bytes: int = 500):
    """Create a random text file for encoding"""
    content = ''.join(random.choices(string.ascii_letters + string.digits + ' \n', k=size_bytes))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(content)
    return output_path

#function that test if the current atributes if they decode succssesfuly
def one_sim(encoding, params):
    # Run the encoding and decoding simulation for the given encoding method
    # This is a placeholder for the actual simulation logic
    # Replace this with the actual implementation of your encoding and decoding process
    try:
        sim = Simulation(params)
        run = sim.run(paralel=False)
        for key, data in run.items():
            if run[key].get('status') == 'SUCCESS':
                print(f"Simulation for {encoding} succeeded.")
                return True
            elif run[key].get('status') == 'FAILURE':
                raise RuntimeError(f"Simulation for {encoding} failed with error: {run[key].get('error')}")
    except Exception as e:
        print(f"Error during simulation for {encoding}: {e}")
        return False


def default_values_of_encoding(encoding):
    if encoding == 'yinyang':
        return {
            'name': 'yy_default',
            'sequence_length': 200,
            'encoding_method': 'yinyang',
            'max_homopolymer': 3,
            'max_content': 0.6,
        }
    if encoding == 'goldman':
        return {
            'name': 'goldman_default',
            'sequence_length': 200,
            'encoding_method': 'goldman',
        }
    if encoding == 'church':
        return {
            'name': 'church_default',
            'sequence_length': 200,
            'encoding_method': 'church',
            'rs_num': 5,
            'max_homopolymer': 3,
            'add_redundancy': False,
            'add_primer': False,
        }
    if encoding == 'gcplus':
        return {
            'name': 'gcplus_default',
            'sequence_length': 200,
            'encoding_method': 'gcplus',
            'gcplus_k': 168,
            'gcplus_l': 8,
            'gcplus_c1': 8,
            'barcode_length': 0,
            'left_primer': '',
            'right_primer': '',
        }
    if encoding == 'max_density':
        return {
            'name': 'max_density_default',
            'codeword_length': 200,
            'sequence_length': 200,
            'encoding_method': 'max_density',
            'max_homopolymer': 3,
            'dna_barcode_length': 8,
            'codeword_maxlength_positions': 8,
            'outer_error_correction': 'reedsolomon',
            'inner_error_correction': None,
            'reed_solo_percentage': 0.8,
            'percent_of_symbols': 1,
            'ltcode_header': 20,
        }
    if encoding == 'no_homopolymer':
        return {
            'name': 'no_homopolymer_default',
            'codeword_length': 200,
            'sequence_length': 200,
            'encoding_method': 'no_homopolymer',
            'max_homopolymer': 3,
            'dna_barcode_length': 8,
            'codeword_maxlength_positions': 8,
            'outer_error_correction': 'reedsolomon',
            'inner_error_correction': None,
            'reed_solo_percentage': 0.8,
            'percent_of_symbols': 1,
            'ltcode_header': 20,
        }
    if encoding == 'wukong':
        return {
            'name': 'wukong_default',
            'sequence_length': 200,
            'encoding_method': 'wukong',
            'max_homopolymer': 3,
            'min_gc_content': 0.4,
            'max_gc_content': 0.6,
            'rule_num': 1,
            'rs_num': 5,
            'add_redundancy': False,
            'add_primer': False,
        }
    if encoding == 'hedges':
        return {
            'name': 'hedges_default',
            'sequence_length': 200,
            'encoding_method': 'hedges',
            'max_homopolymer': 3,
            'gc_max': 8,
            'gc_window': 12,
            'hedges_coderate': 3,
            'kmer_size_cluster': 180,
        }

#funtion that ecnodes and then simulates nanopores and counts erros
def nanopore_simulation(encoding, params):
    # Run the encoding and decoding simulation for the given encoding method
    # This is a placeholder for the actual simulation logic
    # Replace this with the actual implementation of your encoding and decoding process
    path = create_random_file('./tests/testfiles/random_input.txt', size_bytes=500)
    params.filename = path  # Set the filename in params to the generated random file
    bin = Binarize(params)

    # Handle both file_paths (for compressed) and filename (for default binarization)
    # Prepend ./tests/testfiles/ like testbase does
    if hasattr(params, 'file_paths') and params.file_paths:
        file_paths = ['./tests/testfiles/' + fp for fp in params.file_paths]
    elif hasattr(params, 'filename'):
        file_paths = [params.filename]
    else:
        raise ValueError("params must have either 'file_paths' or 'filename' attribute")
    
    data_obj = Data(file_paths=file_paths)
    binary_code = bin.binarize(data_obj)

    coder = Encode(params)
    data_enc, info = coder.encode(binary_code)

    syn = SimulateSynthesis(params)
    data_syn, info = syn.simulate(data_enc)

    sto = SimulateStorage(params)
    data_sto, info = sto.simulate(data_syn)

    seq = SimulateSequencing(params)
    data_seq, info = seq.simulate(data_sto)
    data_seq = InSilicoDNA(data_seq.data)

    print(data_seq.data)


if __name__ == '__main__':
    encodings = ['goldman', 'church', 'gcplus', 'max_density', 'no_homopolymer', 'wukong', 'yinyang', 'hedges']

    base_params = {
        'filename': 'Bohemian_Rhapsody_Lyrics.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': None,
        'mean': 10,
        'recovery_method': 'debruijn_fixedlength',
        'min_coverage': 1,
        'clustering_method': 'kmere_cluster',
        'storage_conditions': None,
        'kmer_seed': 42,   
    }

    nanopore_ids = [39, 40]



    defoultparams = {}      
    for encoding in encodings:
        parameters=base_params.copy()
        parameters.update(default_values_of_encoding(encoding))
        print(f"Parameters for {encoding}: {parameters}")
        defoultparams[encoding] = Params(**parameters)
        
    sim_resoults = {}
    for encoding in encodings:
        sim_resoults[encoding] = one_sim(encoding, defoultparams[encoding])

    print("Simulation results for different error rates:")

    breakpoint()  # Add a breakpoint here for debugging
    
    for encoding in encodings:
        sim_resoults[encoding] = one_sim(encoding, defoultparams[encoding])

    for encoding, error_data in sim_resoults.items():
        print(f"\nEncoding: {encoding}")
        for error_rate, runs in error_data.items():
            success_count = sum(run['success'] for run in runs)
            print(f"  Error Rate: {error_rate}, Successful Runs: {success_count}/{len(runs)}")

