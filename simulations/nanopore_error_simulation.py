import os
import string
import random
import numpy as np
import matplotlib.pyplot as plt

from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.encode import Encode
from dnabyte.params import Params
from dnabyte.sequence import SimulateSequencing
from dnabyte.store import SimulateStorage
from dnabyte.synthesize import SimulateSynthesis
from simulations.simulation import Simulation
from dnabyte.binarize import Binarize
from dnabyte.data_classes import Data


# ============================================================
# Create a random input file
# ============================================================

def create_random_file(output_path: str, size_bytes: int = 500):
    """Create a random text file for encoding."""
    content = ''.join(
        random.choices(
            string.ascii_letters + string.digits + ' \n',
            k=size_bytes
        )
    )

    directory = os.path.dirname(output_path)
    if directory:
        os.makedirs(directory, exist_ok=True)

    with open(output_path, 'w') as f:
        f.write(content)

    return output_path


# ============================================================
# Check whether encoding/decoding succeeds
# ============================================================

def one_sim(encoding, params):

    try:
        print(f"Running simulation for encoding: {encoding}")

        sim = Simulation([params])
        run = sim.run(paralel=False)

        for key, result in run.items():

            if result.get('status') == 'SUCCESS':
                print(f"Simulation for {encoding} succeeded.")
                return True

            elif result.get('status') == 'FAILURE':
                raise RuntimeError(
                    f"Simulation for {encoding} failed: "
                    f"{result.get('error')}"
                )

    except Exception as e:
        print(f"Error during simulation for {encoding}: {e}")
        return False


# ============================================================
# Default parameters for each encoding
# ============================================================

def default_values_of_encoding(encoding):

    if encoding == 'yinyang':
        return {
            'name': 'yy_default',
            'sequence_length': 200,
            'encoding_method': 'yinyang',
            'max_homopolymer': 3,
            'max_content': 0.6,
        }

    elif encoding == 'goldman':
        return {
            'name': 'goldman_default',
            'sequence_length': 200,
            'encoding_method': 'goldman',
        }

    elif encoding == 'church':
        return {
            'name': 'church_default',
            'sequence_length': 200,
            'encoding_method': 'church',
            'rs_num': 5,
            'max_homopolymer': 3,
            'add_redundancy': False,
            'add_primer': False,
        }

    elif encoding == 'gcplus':
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

    elif encoding == 'max_density':
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

    elif encoding == 'no_homopolymer':
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

    elif encoding == 'wukong':
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

    elif encoding == 'hedges':
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

    else:
        raise ValueError(f"Unknown encoding: {encoding}")


# ============================================================
# Nanopore sequencing simulation
# ============================================================

def nanopore_simulation(encoding, params):

    path = create_random_file(
        './tests/testfiles/random_input.txt',
        size_bytes=500
    )

    params.filename = path

    binarizer = Binarize(params)

    if hasattr(params, 'file_paths') and params.file_paths:
        file_paths = [
            './tests/testfiles/' + fp
            for fp in params.file_paths
        ]

    elif hasattr(params, 'filename'):
        file_paths = [params.filename]

    else:
        raise ValueError(
            "params must have either 'file_paths' or 'filename'"
        )

    # --------------------------------------------------------
    # Binarization
    # --------------------------------------------------------

    data_obj = Data(file_paths=file_paths)
    binary_code = binarizer.binarize(data_obj)

    # --------------------------------------------------------
    # Encoding
    # --------------------------------------------------------

    coder = Encode(params)
    data_enc, info = coder.encode(binary_code)

    # --------------------------------------------------------
    # Synthesis
    # --------------------------------------------------------

    syn = SimulateSynthesis(params)
    data_syn, info = syn.simulate(data_enc)

    # --------------------------------------------------------
    # Storage
    # --------------------------------------------------------

    sto = SimulateStorage(params)
    data_sto, info = sto.simulate(data_syn)

    # --------------------------------------------------------
    # Sequencing
    # --------------------------------------------------------

    seq = SimulateSequencing(params)
    data_seq, info = seq.simulate(data_sto)

    # --------------------------------------------------------
    # Calculate total number of sequenced base pairs
    # --------------------------------------------------------

    total_bp = 0

    for strand in data_seq.data:

        # If the strands are strings
        if isinstance(strand, str):
            total_bp += len(strand)

        # If they are objects with a sequence attribute
        elif hasattr(strand, 'sequence'):
            total_bp += len(strand.sequence)

        # Fallback
        else:
            total_bp += len(strand)

    data_seq = InSilicoDNA(data_seq.data)

    return data_seq, info, total_bp


# ============================================================
# Main program
# ============================================================

if __name__ == '__main__':

    encodings = [
        'goldman',
        'church',
        'gcplus',
        'max_density',
        'no_homopolymer',
        'wukong',
        'yinyang',
        'hedges'
    ]

    base_params = {
        'filename': 'Bohemian_Rhapsody_Lyrics.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': None,
        'mean': 1,
        'recovery_method': 'pass_through',
        'min_coverage': 1,
        'clustering_method': 'pass_through',
        'storage_conditions': None,
        'kmer_seed': 42,
    }

    nanopore_ids = [39, 40]

    defaultparams = {}

    # ========================================================
    # Create Params objects
    # ========================================================

    for encoding in encodings:

        parameters = base_params.copy()
        parameters.update(default_values_of_encoding(encoding))

        print(f"\nParameters for {encoding}:")
        print(parameters)

        defaultparams[encoding] = Params(**parameters)

    # ========================================================
    # Test successful encoding/decoding
    # ========================================================

    simulation_results = {}

    for encoding in encodings:
        simulation_results[encoding] = one_sim(
            encoding,
            defaultparams[encoding]
        )

    print("\nSimulation results:")
    print(simulation_results)
# ============================================================
# Nanopore simulations - BOTH IDs
# ============================================================

nanopore_ids = [39, 40]
num_runs = 50

# Structure:
#
# percent_error_counts_per_error_type[39]["goldman"]["insertions"]
# percent_error_counts_per_error_type[40]["goldman"]["insertions"]
#
percent_error_counts_per_error_type = {
    39: {},
    40: {}
}

# Keep the raw totals too, if you want them later
nanopore_totals = {
    39: {},
    40: {}
}

info = {
    39: {},
    40: {}
}


for nanopore_id in nanopore_ids:

    print("\n========================================")
    print(f"RUNNING NANOPore ID {nanopore_id}")
    print("========================================")

    for encoding in encodings:

        print(f"\nRunning {encoding} with nanopore ID {nanopore_id}")

        # ----------------------------------------------------
        # Set nanopore parameters
        # ----------------------------------------------------

        defaultparams[encoding].sequencing_method = 'mesa'
        defaultparams[encoding].mesa_sequencing_mean = 10
        defaultparams[encoding].mesa_sequencing_id = nanopore_id

        # Recreate Params object
        defaultparams[encoding] = Params(
            **defaultparams[encoding].__dict__
        )

        # ----------------------------------------------------
        # Accumulate results from 5 runs
        # ----------------------------------------------------

        total_insertions = 0
        total_deletions = 0
        total_substitutions = 0
        total_bp = 0

        for j in range(num_runs):

            print(
                f"  Run {j + 1}/{num_runs}"
            )

            _, run_info, run_total_bp = nanopore_simulation(
                encoding,
                defaultparams[encoding]
            )

            # Save info
            info[nanopore_id][encoding] = run_info

            # Add errors
            total_insertions += run_info.get(
                'total_insertions', 0
            )

            total_deletions += run_info.get(
                'total_deletions', 0
            )

            total_substitutions += run_info.get(
                'total_substitutions', 0
            )

            # Add total sequenced bp
            total_bp += run_total_bp

        # ----------------------------------------------------
        # Store raw totals
        # ----------------------------------------------------

        nanopore_totals[nanopore_id][encoding] = {
            'insertions': total_insertions,
            'deletions': total_deletions,
            'substitutions': total_substitutions,
            'total_bp': total_bp
        }

        # ----------------------------------------------------
        # Calculate percentage relative to TOTAL bp
        # ----------------------------------------------------

        if total_bp > 0:

            percent_error_counts_per_error_type[
                nanopore_id
            ][encoding] = {

                'insertions': (
                    total_insertions / total_bp
                ) * 100,

                'deletions': (
                    total_deletions / total_bp
                ) * 100,

                'substitutions': (
                    total_substitutions / total_bp
                ) * 100
            }

        else:

            percent_error_counts_per_error_type[
                nanopore_id
            ][encoding] = {
                'insertions': 0,
                'deletions': 0,
                'substitutions': 0
            }


# ============================================================
# Print results
# ============================================================

for nanopore_id in nanopore_ids:

    print(
        f"\n========================================"
    )

    print(
        f"NANOpore ID {nanopore_id}"
    )

    print(
        f"========================================"
    )

    for encoding in encodings:

        percentages = (
            percent_error_counts_per_error_type[
                nanopore_id
            ][encoding]
        )

        print(f"\n{encoding}")

        print(
            f"  Insertions: "
            f"{percentages['insertions']:.4f}%"
        )

        print(
            f"  Deletions: "
            f"{percentages['deletions']:.4f}%"
        )

        print(
            f"  Substitutions: "
            f"{percentages['substitutions']:.4f}%"
        )

# ============================================================
# Plot: Nanopore ID 39
# ============================================================

x = np.arange(len(encodings))
width = 0.25

insertions_39 = [
    percent_error_counts_per_error_type[39][e]['insertions']
    for e in encodings
]

deletions_39 = [
    percent_error_counts_per_error_type[39][e]['deletions']
    for e in encodings
]

substitutions_39 = [
    percent_error_counts_per_error_type[39][e]['substitutions']
    for e in encodings
]


fig, ax = plt.subplots(figsize=(14, 7))

ax.bar(
    x - width,
    insertions_39,
    width,
    label='Insertions'
)

ax.bar(
    x,
    deletions_39,
    width,
    label='Deletions'
)

ax.bar(
    x + width,
    substitutions_39,
    width,
    label='Substitutions'
)

ax.set_xlabel('Encoding Method')

ax.set_ylabel(
    'Errors / Total Sequenced bp (%)'
)

ax.set_title(
    'Nanopore Sequencing Error Rates - ID 39'
)

ax.set_xticks(x)
ax.set_xticklabels(
    encodings,
    rotation=45,
    ha='right'
)

ax.legend()

plt.tight_layout()
plt.show()


# ============================================================
# Plot: Nanopore ID 40
# ============================================================

insertions_40 = [
    percent_error_counts_per_error_type[40][e]['insertions']
    for e in encodings
]

deletions_40 = [
    percent_error_counts_per_error_type[40][e]['deletions']
    for e in encodings
]

substitutions_40 = [
    percent_error_counts_per_error_type[40][e]['substitutions']
    for e in encodings
]


fig, ax = plt.subplots(figsize=(14, 7))

ax.bar(
    x - width,
    insertions_40,
    width,
    label='Insertions'
)

ax.bar(
    x,
    deletions_40,
    width,
    label='Deletions'
)

ax.bar(
    x + width,
    substitutions_40,
    width,
    label='Substitutions'
)

ax.set_xlabel('Encoding Method')

ax.set_ylabel(
    'Errors / Total Sequenced bp (%)'
)

ax.set_title(
    'Nanopore Sequencing Error Rates - ID 40'
)

ax.set_xticks(x)
ax.set_xticklabels(
    encodings,
    rotation=45,
    ha='right'
)

ax.legend()

plt.tight_layout()
plt.show()
