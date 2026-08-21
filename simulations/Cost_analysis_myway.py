
import os

from dnabyte.params import Params
from dnabyte.binarize import Binarize
from simulations.simulation import Simulation
from dnabyte.data_classes.base import Data
from dnabyte.encode import Encode


#Function to tun all encoding with out error correction an finding the highest strand/bp count; this is baseline
def encode_and_count(encodings, defultparams):
    encodedbp = {}
    encsuccess = {}
    for encoding in encodings:
        #sanity check for encoding method
        sim = Simulation([defultparams[encoding]])
        results = sim.run()
        
        # Check if any simulation succeeded
        encsuccess[encoding] = any(
            sim_result.get('status') == 'SUCCESS' 
            for sim_result in results.values()
        )

        bin = Binarize(defultparams[encoding])
        print(f"Encoding: {encoding}, Params: {defultparams[encoding]}")
        print(defultparams[encoding].filename)
        
        # Handle both file_paths (for compressed) and filename (for default binarization)
        # Prepend ./tests/testfiles/ like testbase does
        if hasattr(defultparams[encoding], 'file_paths') and defultparams[encoding].file_paths:
            file_paths = ['./tests/testfiles/' + fp for fp in defultparams[encoding].file_paths]
        elif hasattr(defultparams[encoding], 'filename'):
            file_paths = ['./tests/testfiles/' + defultparams[encoding].filename]
        else:
            raise ValueError("params must have either 'file_paths' or 'filename' attribute")
        
        data_obj = Data(file_paths=file_paths)
        binary_code = bin.binarize(data_obj)
        coder = Encode(defultparams[encoding])
        data_enc, info = coder.encode(binary_code)
        encodedbp[encoding] = sum(len(i) for i in data_enc.data)  # Sum lengths of all encoded sequences

    print(f"Encoded bp counts: {encodedbp}")
    print(f"Encoding success status: {encsuccess}")

    if all(value is True for value in encsuccess.values()):
        return encodedbp
    else:
        failed_encodings = [encoding for encoding, success in encsuccess.items() if not success]
        raise RuntimeError(f"Baseline simulation failed for the following encodings: {', '.join(failed_encodings)}")

#Function with the found amount of bp that increses the enecoding parameters of all other encodings to have the same amount of bp and then run all encodings with error correction and find the highest strand/bp count;

# main function run for the found highest strand/bp count for diffrent error rates and find the highest strand/bp count for each encoding with error correction and without error correction;


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
            'rs_num': 0,
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
            'gcplus_c1': 0,
            'barcode_length': 0,
            'left_primer': '',
            'right_primer': '',
        }
    if encoding == 'max_density':
        return {
            'name': 'max_density_default',
            'sequence_length': 200,
            'encoding_method': 'max_density',
            'max_homopolymer': 3,
            'dna_barcode_length': 8,
            'codeword_maxlength_positions': 8
        }
    if encoding == 'no_homopolymer':
        return {
            'name': 'no_homopolymer_default',
            'sequence_length': 200,
            'encoding_method': 'no_homopolymer',
            'max_homopolymer': 3,
            'dna_barcode_length': 8,
            'codeword_maxlength_positions': 8
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
            'rs_num': 0,
            'add_redundancy': False,
            'add_primer': False,
        }
    if encoding == 'hedges':
        return {
            'name': 'hedges_default',
            'sequence_length': 200,
            'encoding_method': 'hedges',
            'max_homopolymer': 3,
            'gc_max': 12,
            'gc_window': 8,
            'hedges_coderate': 3,
            'kmer_size_cluster': 180,
        }


if __name__ == '__main__':
    encodings = ['goldman', 'church', 'gcplus', 'max_density', 'no_homopolymer', 'wukong', 'yinyang', 'hedges']

    base_params = {
        'filename': 'Bohemian_Rhapsody_Lyrics.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'kmere',
        'kmer_k': 1,
        'recovery_method': 'debruijn',
        'min_coverage': 1,
        'clustering_method': 'kmere_cluster',
        'storage_conditions': None,
        'kmer_seed': 42,
    }

    no_error_params = {
        'mean': 1,
        'std_dev': 0,
        'p_ins': 0,
        'p_del': 0,
        'p_sub': 0,
    }
    
    error_rates = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.12, 0.16, 0.20]

    defoultparams = {}
    for encoding in encodings:
        parameters=base_params.copy()
        parameters.update(default_values_of_encoding(encoding))
        defoultparams[encoding] = Params(**parameters)
        

    encode_and_count(encodings, defoultparams)
    
    