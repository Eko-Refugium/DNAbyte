import copy
import random

import matplotlib.pyplot as plt

from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.data_classes.base import Data
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.encode import Encode
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.misc_err import SimulateMiscErrors
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus


# ---------------------------------------------------------------------------
# Keep the same default style as the original script, but make evaluation
# stable by encoding once and then reusing the same encoded strands.
# ---------------------------------------------------------------------------

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
            'codeword_length': 200,
            'sequence_length': 200,
            'encoding_method': 'max_density',
            'max_homopolymer': 3,
            'dna_barcode_length': 8,
            'codeword_maxlength_positions': 8,
            'outer_error_correction': None,
            'inner_error_correction': None,
            'reed_solo_percentage': 1,
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
            'outer_error_correction': None,
            'inner_error_correction': None,
            'reed_solo_percentage': 1,
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
            'gc_max': 8,
            'gc_window': 12,
            'hedges_coderate': 3,
            'kmer_size_cluster': 180,
        }

    raise ValueError(f'Unknown encoding: {encoding}')


def build_default_params(encodings):
    base_params = {
        'filename': 'Bohemian_Rhapsody_Lyrics.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'iid',
        'recovery_method': 'simple',
        'min_coverage': 1,
        'clustering_method': 'kmere_cluster',
        'storage_conditions': None,
        'kmer_seed': 42,
        'mean': 1,
        'std_dev': 0,
        'iid_error_rate': 0.0,
        'iid_substitution_rate': 0.0,
        'iid_insertion_rate': 0.0,
        'iid_deletion_rate': 0.0,
    }

    params = {}
    for encoding in encodings:
        current = base_params.copy()
        current.update(default_values_of_encoding(encoding))
        params[encoding] = Params(**current)
    return params


def binarize_data_for_params(params):
    bin_obj = Binarize(params)

    if hasattr(params, 'file_paths') and params.file_paths:
        file_paths = ['./tests/testfiles/' + fp for fp in params.file_paths]
    elif hasattr(params, 'filename'):
        file_paths = ['./tests/testfiles/' + params.filename]
    else:
        raise ValueError("params must have either 'file_paths' or 'filename' attribute")

    data_obj = Data(file_paths=file_paths)
    return bin_obj.binarize(data_obj)


def encode_once_for_all(encodings, default_params):
    encoded = {}

    for encoding in encodings:
        params = default_params[encoding]
        binary_code = binarize_data_for_params(params)
        coder = Encode(params)
        encoded_data, info = coder.encode(binary_code)

        encoded[encoding] = {
            'params': params,
            'binary_code': binary_code,
            'encoded_data': encoded_data,
            'info': info,
        }

    return encoded


def get_param_snapshot(params, overrides=None):
    data = params.__dict__.copy()
    if overrides:
        data.update(overrides)
    return Params(**data)


def evaluate_single_fixed_encoding(encoding, params_template, binary_code, encoded_data, error_rate, trial_index):
    # Rebuild a fresh Params object for this particular trial so you are not mutating
    # the same config object across runs.
    params = get_param_snapshot(
        params_template,
        {
            'seed': params_template.kmer_seed + trial_index + int(error_rate * 1000),
            'iid_error_rate': error_rate,
            'iid_substitution_rate': error_rate,
            'iid_insertion_rate': 0.0,
            'iid_deletion_rate': 0.0,
        }
    )

    # Use the same encoded object for every trial at that rate.
    current_data = encoded_data
    if not isinstance(current_data, NucleobaseCode):
        current_data = NucleobaseCode(current_data.data)

    if params.synthesis_method is not None:
        current_data, info = SimulateSynthesis(params).simulate(current_data)

    if params.storage_conditions is not None:
        current_data, info = SimulateStorage(params).simulate(current_data)

    if params.error_methods is not None:
        current_data, info = SimulateMiscErrors(params).simulate(current_data)

    if params.sequencing_method is not None:
        current_data, info = SimulateSequencing(params).simulate(current_data)
        current_data = InSilicoDNA(current_data.data)

    if hasattr(params, 'clustering_method') and hasattr(params, 'recovery_method') and params.clustering_method and params.recovery_method:
        cluster_obj = Cluster(params)
        current_data, info = cluster_obj.cluster(current_data)
        consensus_obj = Consensus(params)
        current_data, info = consensus_obj.call(current_data)
    else:
        coder = Encode(params)
        current_data, info = coder.process(current_data)

    coder = Encode(params)
    data_dec, valid, info = coder.decode(current_data)

    try:
        comparison, res = data_dec.compare(data_dec, binary_code)
        success = (comparison != 'ERROR') and valid
    except Exception:
        success = bool(valid)

    return success


def simulate_cost_analysis_once_encoded(encodings, default_params, encoded_once, error_rates, trials=25):
    results = {}

    for encoding in encodings:
        results[encoding] = {}
        for error_rate in error_rates:
            run_successes = []
            for trial in range(trials):
                success = evaluate_single_fixed_encoding(
                    encoding,
                    encoded_once[encoding]['params'],
                    encoded_once[encoding]['binary_code'],
                    encoded_once[encoding]['encoded_data'],
                    error_rate,
                    trial,
                )
                run_successes.append(success)
            results[encoding][error_rate] = run_successes

    return results


def print_summary(results):
    for encoding, error_data in results.items():
        print(f"\nEncoding: {encoding}")
        for error_rate, runs in sorted(error_data.items()):
            success_count = sum(runs)
            print(f"  Error Rate: {error_rate}, Successful Runs: {success_count}/{len(runs)}")


if __name__ == '__main__':
    encodings = ['yinyang', 'hedges', 'wukong', 'no_homopolymer', 'church', 'max_density', 'goldman']
    default_params = build_default_params(encodings)
    encoded_once = encode_once_for_all(encodings, default_params)

    error_rates = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2]

    sim_results = simulate_cost_analysis_once_encoded(
        encodings,
        default_params,
        encoded_once,
        error_rates,
        trials=20,
    )

    print_summary(sim_results)

    for encoding, error_data in sim_results.items():
        error_rates_sorted = sorted(error_data.keys())
        success_rates = [sum(error_data[er]) / len(error_data[er]) for er in error_rates_sorted]
        plt.plot(error_rates_sorted, success_rates, label=encoding)

    plt.xlabel('Error Rate')
    plt.ylabel('Success Rate')
    plt.title('Simulation Results for Different Error Rates with a single shared encoding per method')
    plt.legend()
    plt.tight_layout()
    plt.show()
