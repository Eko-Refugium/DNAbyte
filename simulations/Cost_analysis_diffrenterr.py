import copy
import csv
import json
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np

from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.data_classes import Data
from dnabyte.encode import Encode
from dnabyte.sequence import SimulateSequencing
from dnabyte.store import SimulateStorage
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus


# ================================================================
# E2 — SEQUENCING ERROR COMPOSITION
# ================================================================
#
# Structure follows E1:
#   encoding -> configuration -> error profile -> error rate -> repeats
#
# Fixed physical redundancy:
#   mean = 100
#   mmseq_min = 20
#
# E2 changes the COMPOSITION of a fixed total sequencing error rate.
#
# Error profiles:
#   substitution
#   insertion
#   deletion
#   substitution + insertion
#   substitution + deletion
#   insertion + deletion
#   balanced
#
# Outputs:
#   * raw JSON
#   * summary CSV
#   * recovery-vs-error plots
#   * heatmaps
#   * scatter plot
#   * bar plot at 1% total error
#   * DNA-cost-vs-tolerance plots
#
# ================================================================


# ================================================================
# EXPERIMENT SETTINGS
# ================================================================

ENCODINGS = [
    "goldman",
    "church",
    "no_homopolymer",
    "max_density",
    "gcplus",
    "hedges",
    "wukong",
    "yinyang",
]


REPEATS = 100

FIXED_MEAN = 1

MMSEQ_MIN = 20


ERROR_RATES = [
    0.0,
    0.0001,
    0.0002,
    0.0003,
    0.0004,
    0.0005,
    0.0006,
    0.0007,
    0.0008,
    0.0009,
    0.001,
    0.002,
    0.003,
    0.004,
    0.005,
    0.006,
    0.007,
    0.008,
    0.009,
    0.01,
    0.02,
    0.03,
    0.04,
    0.05,
    0.06,
    0.07,
    0.08,
    0.09,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
]


ERROR_PROFILES = {
    "Substitution": (1.0, 0.0, 0.0),
    "Insertion": (0.0, 1.0, 0.0),
    "Deletion": (0.0, 0.0, 1.0),
    "Sub + Ins": (0.5, 0.5, 0.0),
    "Sub + Del": (0.5, 0.0, 0.5),
    "Ins + Del": (0.0, 0.5, 0.5),
    "Balanced": (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
}


# Keep this True to use the same configuration sweep as E1.
# Set False for a much smaller first test using one baseline
# configuration per encoding.
USE_E1_CONFIGURATIONS = True


# ================================================================
# INPUT / OUTPUT
# ================================================================

INPUT_FILENAME = "testfilesimsall.txt"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

OUTPUT_DIR = os.path.join(
    SCRIPT_DIR,
    "simulation_results_E2"
)

os.makedirs(OUTPUT_DIR, exist_ok=True)

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


# ================================================================
# E1 ENCODING DEFAULTS
# ================================================================


def default_values_of_encoding(encoding):

    if encoding == "yinyang":
        return {
            "name": "yy_default",
            "sequence_length": 200,
            "encoding_method": "yinyang",
            "max_homopolymer": 3,
            "max_content": 0.6,
        }

    if encoding == "goldman":
        return {
            "name": "goldman_default",
            "sequence_length": 200,
            "encoding_method": "goldman",
            "kmer_size_cluster": 180,
        }

    if encoding == "church":
        return {
            "name": "church_default",
            "sequence_length": 200,
            "encoding_method": "church",
            "rs_num": 0,
            "max_homopolymer": 3,
            "add_redundancy": False,
            "add_primer": False,
        }

    if encoding == "gcplus":
        return {
            "name": "gcplus_default",
            "sequence_length": 200,
            "encoding_method": "gcplus",
            "gcplus_k": 168,
            "gcplus_l": 8,
            "gcplus_c1": 2,
            "barcode_length": 0,
            "left_primer": "",
            "right_primer": "",
            "kmer_size_debruijn": 90,
            "kmer_size_cluster": 20,
            "kmer_threshold": 0.9,
        }

    if encoding == "max_density":
        return {
            "name": "max_density_default",
            "codeword_length": 200,
            "sequence_length": 200,
            "encoding_method": "max_density",
            "max_homopolymer": 3,
            "dna_barcode_length": 8,
            "codeword_maxlength_positions": 8,
            "outer_error_correction": None,
            "inner_error_correction": None,
            "reed_solo_percentage": 1,
            "percent_of_symbols": 1,
            "ltcode_header": 20,
        }

    if encoding == "no_homopolymer":
        return {
            "name": "no_homopolymer_default",
            "codeword_length": 200,
            "sequence_length": 200,
            "encoding_method": "no_homopolymer",
            "max_homopolymer": 3,
            "dna_barcode_length": 8,
            "codeword_maxlength_positions": 8,
            "outer_error_correction": None,
            "inner_error_correction": None,
            "reed_solo_percentage": 1,
            "percent_of_symbols": 1,
            "ltcode_header": 20,
        }

    if encoding == "wukong":
        return {
            "name": "wukong_default",
            "sequence_length": 200,
            "encoding_method": "wukong",
            "max_homopolymer": 3,
            "min_gc_content": 0.4,
            "max_gc_content": 0.6,
            "rule_num": 1,
            "rs_num": 0,
            "add_redundancy": False,
            "add_primer": False,
        }

    if encoding == "hedges":
        return {
            "name": "hedges_default",
            "sequence_length": 200,
            "encoding_method": "hedges",
            "max_homopolymer": 3,
            "gc_max": 8,
            "gc_window": 12,
            "hedges_coderate": 3,
            "kmer_size_cluster": 180,
        }

    raise ValueError(f"Unknown encoding: {encoding}")


# ================================================================
# E1 PARAMETER CONFIGURATIONS
# ================================================================


def get_parameter_configurations(encoding):

    if encoding == "church":
        return [
            {
                "add_redundancy": add_redundancy,
                "rs_num": rs_num,
            }
            for add_redundancy in [False, True]
            for rs_num in range(0, 6)
        ]

    if encoding == "wukong":
        return [
            {
                "add_redundancy": add_redundancy,
                "rs_num": rs_num,
            }
            for add_redundancy in [False, True]
            for rs_num in range(0, 6)
        ]

    if encoding == "goldman" or encoding == "yinyang":
        return [
            {"mean": 1}
        ]

    if encoding == "gcplus":
        return [
            {"gcplus_c1": c1}
            for c1 in [1, 2, 3, 4, 5]
        ]

    if encoding == "hedges":
        return [
            {"hedges_coderate": rate}
            for rate in [1, 2, 3, 4, 5, 6]
        ]

    if encoding in ["max_density", "no_homopolymer"]:
        return [
            {"reed_solo_percentage": percentage}
            for percentage in [
                0.50,
                0.60,
                0.70,
                0.80,
                0.90,
            ]
        ]

    raise ValueError(
        f"No parameter sweep defined for {encoding}"
    )


# ================================================================
# BASE PARAMETER CREATION
# ================================================================


def build_base_params():

    base_params = {
        "filename": INPUT_FILENAME,
        "binarization_method": "default",
        "synthesis_method": "nosynthpoly",
        "sequencing_method": "iid",
        "recovery_method": "pass_through",
        "min_coverage": 1,
        "clustering_method": "pass_through",
        "storage_conditions": None,
        "kmer_seed": 42,
        "similarity_threshold": 1,
        "gap_penalty": 1.0,

        # Fixed physical redundancy for E2.
        "mean": FIXED_MEAN,
        "std_dev": 0,

        # E1 explicitly controls the three IID components.
        "iid_error_rate": 0.0,
        "iid_substitution_rate": 0.0,
        "iid_insertion_rate": 0.0,
        "iid_deletion_rate": 0.0,
    }

    return base_params


# ================================================================
# CONFIGURATION LIST
# ================================================================


def get_configurations(encoding):

    if USE_E1_CONFIGURATIONS:
        return get_parameter_configurations(encoding)

    return [{}]


# ================================================================
# CONFIGURATION LABEL
# ================================================================


def configuration_name(config):

    if not config:
        return "baseline"

    return "_".join(
        f"{key}={value}"
        for key, value in config.items()
    )


# ================================================================
# ERROR PROFILE LABEL
# ================================================================


def error_profile_rates(total_error_rate, profile):

    sub_fraction, ins_fraction, del_fraction = profile

    return {
        "substitution": total_error_rate * sub_fraction,
        "insertion": total_error_rate * ins_fraction,
        "deletion": total_error_rate * del_fraction,
    }


# ================================================================
# SET E2 SEQUENCING ERROR PARAMETERS
# ================================================================


def set_error_parameters(params, total_error_rate, profile):

    rates = error_profile_rates(
        total_error_rate,
        profile
    )

    # Follow E1's working behavior: use the component rates directly.
    params.iid_substitution_rate = rates["substitution"]
    params.iid_insertion_rate = rates["insertion"]
    params.iid_deletion_rate = rates["deletion"]

    # E1 does not use iid_error_rate for this sweep.
    params.iid_error_rate = 0.0


# ================================================================
# PREPARE PARAMETER OBJECT
# ================================================================


def make_params(encoding, config):

    parameters = build_base_params()
    parameters.update(default_values_of_encoding(encoding))

    # Apply configuration from E1.
    parameters.update(config)

    params = Params(**parameters)

    # ------------------------------------------------------------
    # E2 FIX: physical redundancy is always 100.
    # ------------------------------------------------------------
    #
    # Goldman's E1 configuration contains mean=1 because E1 was
    # sweeping mean. E2 is specifically a sequencing-error
    # composition experiment, so all encodings use mean=100.
    # ------------------------------------------------------------

    params.mean = FIXED_MEAN
    params.std_dev = 0

    # ------------------------------------------------------------
    # Parameters used by the normal recovery pipeline.
    # ------------------------------------------------------------

    params.sequencing_method = "iid"
    params.storage_conditions = None
    params.synthesis_method = "nosynthpoly"
    params.mmseq_min = MMSEQ_MIN
    params.clustering_method = "mmseq"
    params.recovery_method = "debruijn"

    # These are the settings used for the E2 DNA dataset.
    if hasattr(params, "codeword_length"):
        params.codeword_length = 200
    else:
        params.codeword_length = 200

    return params


# ================================================================
# PREPARE ENCODING-SPECIFIC RECOVERY SETTINGS
# ================================================================


def prepare_special_encoding(params, encoding):

    if encoding in [
        "max_density",
        "no_homopolymer",
    ]:
        params.outer_error_correction = "reedsolomon"
        params.clustering_method = None
        params.recovery_method = None

    else:
        # Normal clustering/recovery path.
        params.clustering_method = "pass_through"
        params.recovery_method = "pass_through"


# ================================================================
# BUILD INPUT DATA + ENCODE ONCE PER CONFIGURATION
# ================================================================


def encode_dataset(params):

    binarizer = Binarize(params)

    if (
        hasattr(params, "file_paths")
        and params.file_paths
    ):
        file_paths = [
            "./tests/testfiles/" + fp
            for fp in params.file_paths
        ]

    elif hasattr(params, "filename"):
        file_paths = [
            "./tests/testfiles/" + params.filename
        ]

    else:
        raise ValueError(
            "params must have either 'file_paths' "
            "or 'filename' attribute"
        )

    data_obj = Data(
        file_paths=file_paths
    )

    binary_code = binarizer.binarize(
        data_obj
    )

    coder = Encode(params)

    data_enc, info = coder.encode(
        binary_code
    )

    if hasattr(data_enc, "data"):
        dna_sequences = data_enc.data
    else:
        dna_sequences = data_enc

    dna_length = sum(
        len(str(sequence))
        for sequence in dna_sequences
    )

    information_bits = len(
        binary_code.data
    )

    dna_cost_nt_per_bit = (
        dna_length / information_bits
        if information_bits > 0
        else 0
    )

    return (
        binary_code,
        data_enc,
        coder,
        information_bits,
        dna_length,
        dna_cost_nt_per_bit,
    )


# ================================================================
# RUN ONE RECOVERY TRIAL
# ================================================================


def run_one_trial(
    params,
    data_enc,
    coder,
    binary_code,
    encoding,
):

    success = False
    comparison_result = None
    error_message = None

    try:

        # --------------------------------------------------------
        # SYNTHESIS
        # --------------------------------------------------------

        syn = SimulateSynthesis(params)

        data_syn, info = syn.simulate(
            data_enc
        )

        # --------------------------------------------------------
        # STORAGE
        # --------------------------------------------------------

        sto = SimulateStorage(params)

        data_sto, info = sto.simulate(
            data_syn
        )

        # --------------------------------------------------------
        # SEQUENCING
        # --------------------------------------------------------

        seq = SimulateSequencing(params)

        data_seq, info = seq.simulate(
            data_sto
        )

        # --------------------------------------------------------
        # RECOVERY / CLUSTERING
        # --------------------------------------------------------

        if (
            params.clustering_method is not None
            and params.recovery_method is not None
        ):

            cluster_obj = Cluster(params)

            data_cluster, info = (
                cluster_obj.cluster(data_seq)
            )

            consensus_obj = Consensus(params)

            data_cor, info = (
                consensus_obj.call(data_cluster)
            )

        else:

            # max_density and no_homopolymer use their own
            # encoding-specific recovery/processing path.

            data_cor, info = coder.process(
                data_seq
            )

        # --------------------------------------------------------
        # DECODE
        # --------------------------------------------------------

        data_dec, valid, info = coder.decode(
            data_cor
        )

        if not valid:
            raise ValueError(
                "decoded data invalid"
            )

        # --------------------------------------------------------
        # COMPARE WITH ORIGINAL
        # --------------------------------------------------------

        comparison_result, result_data = (
            data_dec.compare(
                data_dec,
                binary_code
            )
        )

        success = (
            comparison_result == "SUCCESS"
        )

        return {
            "success": success,
            "comparison": comparison_result,
            "results": result_data,
            "error": None,
        }

    except Exception as exc:

        error_message = (
            f"{type(exc).__name__}: {exc}"
        )

        print(
            f"    {encoding} trial failed: "
            f"{error_message}"
        )

        return {
            "success": False,
            "comparison": None,
            "results": None,
            "error": error_message,
        }


# ================================================================
# MAIN E2 SIMULATION
# ================================================================


def simulate_e2():

    base_params = build_base_params()

    all_results = {}

    total_configurations = sum(
        len(get_configurations(encoding))
        for encoding in ENCODINGS
    )

    configuration_counter = 0

    print()
    print("=" * 90)
    print("E2 — SEQUENCING ERROR COMPOSITION")
    print("=" * 90)
    print()
    print(f"Encodings:             {len(ENCODINGS)}")
    print(f"Configurations:        {total_configurations}")
    print(f"Error profiles:        {len(ERROR_PROFILES)}")
    print(f"Error rates:           {len(ERROR_RATES)}")
    print(f"Repeats:               {REPEATS}")
    print(f"Fixed mean:            {FIXED_MEAN}")
    print(f"MMseq minimum:         {MMSEQ_MIN}")
    print()

    for encoding in ENCODINGS:

        print()
        print("#" * 90)
        print(f"ENCODING: {encoding}")
        print("#" * 90)

        all_results[encoding] = {}

        configurations = get_configurations(
            encoding
        )

        for config_index, config in enumerate(
            configurations,
            start=1
        ):

            configuration_counter += 1

            config_name = configuration_name(
                config
            )

            print()
            print("-" * 90)
            print(
                f"Configuration "
                f"{configuration_counter}/"
                f"{total_configurations}: "
                f"{encoding} / {config_name}"
            )
            print("-" * 90)

            params = make_params(
                encoding,
                config
            )

            prepare_special_encoding(
                params,
                encoding
            )

            # ----------------------------------------------------
            # Encode once for this configuration.
            # Sequencing error composition does not change the
            # encoded DNA or its storage cost.
            # ----------------------------------------------------

            try:

                (
                    binary_code,
                    data_enc,
                    coder,
                    information_bits,
                    dna_length,
                    dna_cost_nt_per_bit,
                ) = encode_dataset(params)

            except Exception as exc:

                print(
                    f"ENCODING FAILED: "
                    f"{type(exc).__name__}: {exc}"
                )

                all_results[encoding][config_name] = {
                    "configuration": copy.deepcopy(config),
                    "information_bits": None,
                    "dna_length": None,
                    "dna_cost_nt_per_bit": None,
                    "profiles": {},
                    "encoding_error": (
                        f"{type(exc).__name__}: {exc}"
                    ),
                }

                continue

            print(
                f"DNA length: {dna_length} nt"
            )

            print(
                f"Information bits: {information_bits}"
            )

            print(
                f"DNA cost: {dna_cost_nt_per_bit:.4f} nt/bit"
            )

            config_results = {
                "configuration": copy.deepcopy(config),
                "information_bits": information_bits,
                "dna_length": dna_length,
                "dna_cost_nt_per_bit": dna_cost_nt_per_bit,
                "profiles": {},
            }

            # ----------------------------------------------------
            # Error composition profiles
            # ----------------------------------------------------

            for profile_name, profile_fractions in ERROR_PROFILES.items():

                print()
                print(
                    f"  ERROR PROFILE: {profile_name}"
                )

                profile_results = {}

                zero_reached = False

                for error_rate in ERROR_RATES:

                    # ------------------------------------------------
                    # Requested early stop:
                    # once 0/100 occurs for this profile/configuration,
                    # all higher error rates are recorded as 0/100.
                    # ------------------------------------------------

                    if zero_reached:

                        profile_results[error_rate] = {
                            "total_error_rate": error_rate,
                            "substitution_rate": (
                                error_rate
                                * profile_fractions[0]
                            ),
                            "insertion_rate": (
                                error_rate
                                * profile_fractions[1]
                            ),
                            "deletion_rate": (
                                error_rate
                                * profile_fractions[2]
                            ),
                            "successful_runs": 0,
                            "total_runs": REPEATS,
                            "recovery_rate": 0.0,
                            "skipped_after_zero": True,
                            "runs": [],
                        }

                        print(
                            f"    {error_rate:<8} "
                            f"0/{REPEATS} "
                            f"(skipped after previous zero)"
                        )

                        continue

                    component_rates = error_profile_rates(
                        error_rate,
                        profile_fractions
                    )

                    # ------------------------------------------------
                    # Set sequencing error composition.
                    # ------------------------------------------------

                    set_error_parameters(
                        params,
                        error_rate,
                        profile_fractions
                    )

                    # ------------------------------------------------
                    # If mean > 1, follow E1's similarity adjustment.
                    # ------------------------------------------------

                    if getattr(params, "mean", 1) > 1:
                        params.similarity_threshold = 0.75

                    expected_errors = (
                        dna_length * error_rate
                    )

                    runs = []

                    successful_runs = 0

                    print(
                        f"    error={error_rate:<8} "
                        f"sub={component_rates['substitution']:<10g} "
                        f"ins={component_rates['insertion']:<10g} "
                        f"del={component_rates['deletion']:<10g}"
                    )

                    # ------------------------------------------------
                    # 100 repeats
                    # ------------------------------------------------

                    for repeat in range(REPEATS):

                        trial = run_one_trial(
                            params,
                            data_enc,
                            coder,
                            binary_code,
                            encoding,
                        )

                        if trial["success"]:
                            successful_runs += 1

                        runs.append({
                            "run": repeat + 1,
                            "success": trial["success"],
                            "comparison": trial["comparison"],
                            "results": trial["results"],
                            "error": trial["error"],
                            "information_bits": information_bits,
                            "dna_length": dna_length,
                            "dna_cost_nt_per_bit": dna_cost_nt_per_bit,
                            "expected_errors": expected_errors,
                            "total_error_rate": error_rate,
                            "substitution_rate": component_rates["substitution"],
                            "insertion_rate": component_rates["insertion"],
                            "deletion_rate": component_rates["deletion"],
                            "encoding": encoding,
                            "configuration": copy.deepcopy(config),
                            "error_profile": profile_name,
                        })

                    total_runs = len(runs)

                    recovery_rate = (
                        successful_runs / total_runs
                        if total_runs > 0
                        else 0.0
                    )

                    profile_results[error_rate] = {
                        "total_error_rate": error_rate,
                        "substitution_rate": component_rates["substitution"],
                        "insertion_rate": component_rates["insertion"],
                        "deletion_rate": component_rates["deletion"],
                        "successful_runs": successful_runs,
                        "total_runs": total_runs,
                        "recovery_rate": recovery_rate,
                        "skipped_after_zero": False,
                        "runs": runs,
                    }

                    print(
                        f"      -> "
                        f"{successful_runs}/{total_runs} "
                        f"({recovery_rate:.1%})"
                    )

                    # ------------------------------------------------
                    # Early stopping requested by the user.
                    # ------------------------------------------------

                    if successful_runs == 0:

                        zero_reached = True

                        print(
                            "      -> 0/100 reached; "
                            "higher error rates will be recorded as 0/100."
                        )

                config_results["profiles"][
                    profile_name
                ] = profile_results

            all_results[encoding][config_name] = config_results

    return all_results


# ================================================================
# SUMMARY HELPERS
# ================================================================


def get_profile_recovery(
    results,
    encoding,
    config_name,
    profile_name,
    error_rate,
):

    return results[
        encoding
    ][
        config_name
    ][
        "profiles"
    ][
        profile_name
    ][
        error_rate
    ][
        "recovery_rate"
    ]



def first_configuration_name(results, encoding):
    return next(iter(results[encoding]))


# ================================================================
# SAVE RAW JSON
# ================================================================


def save_json(results):

    output_file = os.path.join(
        OUTPUT_DIR,
        f"E2_raw_results_{TIMESTAMP}.json"
    )

    payload = {
        "experiment": "E2_sequencing_error_composition",
        "encodings": ENCODINGS,
        "repeats": REPEATS,
        "fixed_mean": FIXED_MEAN,
        "mmseq_min": MMSEQ_MIN,
        "error_rates": ERROR_RATES,
        "error_profiles": ERROR_PROFILES,
        "use_e1_configurations": USE_E1_CONFIGURATIONS,
        "results": results,
    }

    with open(output_file, "w") as f:
        json.dump(
            payload,
            f,
            indent=4,
            default=str,
        )

    print(f"Saved raw results: {output_file}")

    return output_file


# ================================================================
# SAVE SUMMARY CSV
# ================================================================


def save_summary_csv(results):

    output_file = os.path.join(
        OUTPUT_DIR,
        f"E2_summary_{TIMESTAMP}.csv"
    )

    with open(
        output_file,
        "w",
        newline="",
    ) as f:

        writer = csv.writer(f)

        writer.writerow([
            "encoding",
            "configuration",
            "dna_cost_nt_per_bit",
            "error_profile",
            "total_error_rate",
            "substitution_rate",
            "insertion_rate",
            "deletion_rate",
            "successful_runs",
            "total_runs",
            "recovery_rate",
            "skipped_after_zero",
        ])

        for encoding, encoding_data in results.items():

            for config_name, config_data in encoding_data.items():

                if config_data.get("encoding_error"):
                    continue

                dna_cost = config_data[
                    "dna_cost_nt_per_bit"
                ]

                for profile_name, profile_data in config_data[
                    "profiles"
                ].items():

                    for error_rate, point in profile_data.items():

                        writer.writerow([
                            encoding,
                            config_name,
                            dna_cost,
                            profile_name,
                            error_rate,
                            point["substitution_rate"],
                            point["insertion_rate"],
                            point["deletion_rate"],
                            point["successful_runs"],
                            point["total_runs"],
                            point["recovery_rate"],
                            point["skipped_after_zero"],
                        ])

    print(f"Saved summary CSV: {output_file}")

    return output_file


# ================================================================
# SAVE SIMPLE PRINTABLE SUMMARY
# ================================================================


def print_summary(results):

    print()
    print("=" * 90)
    print("E2 SUMMARY")
    print("=" * 90)

    for encoding in ENCODINGS:

        print()
        print(f"ENCODING: {encoding}")

        for config_name, config_data in results[encoding].items():

            print(
                f"  Configuration: {config_name}"
            )

            if config_data.get("encoding_error"):
                print(
                    f"    ERROR: {config_data['encoding_error']}"
                )
                continue

            for profile_name, profile_data in config_data[
                "profiles"
            ].items():

                # Print only the last simulated point for compactness.
                simulated_points = [
                    (rate, point)
                    for rate, point in profile_data.items()
                    if not point["skipped_after_zero"]
                ]

                if not simulated_points:
                    continue

                rate, point = simulated_points[-1]

                print(
                    f"    {profile_name:<15} "
                    f"last_tested={rate:<8g} "
                    f"recovery={point['recovery_rate']:.1%}"
                )


# ================================================================
# PLOT 1 — RECOVERY VS ERROR RATE
# ================================================================
#
# One figure per encoding, using the first E1 configuration as the
# clean reference configuration. Raw results still contain every
# configuration.
#
# This keeps the main paper figure readable rather than producing
# 46 x 7 overlapping lines.
# ================================================================


def plot_reference_recovery_curves(results):

    plot_dir = os.path.join(
        OUTPUT_DIR,
        "plots"
    )

    os.makedirs(plot_dir, exist_ok=True)

    for encoding in ENCODINGS:

        config_name = first_configuration_name(
            results,
            encoding
        )

        config_data = results[
            encoding
        ][
            config_name
        ]

        if config_data.get("encoding_error"):
            continue

        plt.figure(figsize=(12, 7))

        for profile_name in ERROR_PROFILES:

            profile_data = config_data[
                "profiles"
            ][
                profile_name
            ]

            x = list(profile_data.keys())

            y = [
                point["recovery_rate"]
                for point in profile_data.values()
            ]

            plt.plot(
                x,
                y,
                "o-",
                linewidth=2,
                markersize=4,
                label=profile_name,
            )

        plt.xlabel(
            "Total sequencing error rate"
        )

        plt.ylabel(
            "Recovery rate"
        )

        plt.title(
            f"E2: {encoding} — recovery vs error composition\n"
            f"Reference configuration: {config_name}, mean={FIXED_MEAN}"
        )

        plt.xscale("log")
        plt.ylim(0, 1.05)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()

        png_file = os.path.join(
            plot_dir,
            f"E2_{encoding}_recovery_curves.png"
        )

        pdf_file = os.path.join(
            plot_dir,
            f"E2_{encoding}_recovery_curves.pdf"
        )

        plt.savefig(
            png_file,
            dpi=300,
            bbox_inches="tight"
        )

        plt.savefig(
            pdf_file,
            bbox_inches="tight"
        )

        plt.close()

        print(f"Saved: {png_file}")
        print(f"Saved: {pdf_file}")


# ================================================================
# PLOT 2 — HEATMAPS
# ================================================================


def plot_reference_heatmaps(results):

    plot_dir = os.path.join(
        OUTPUT_DIR,
        "plots"
    )

    os.makedirs(plot_dir, exist_ok=True)

    for encoding in ENCODINGS:

        config_name = first_configuration_name(
            results,
            encoding
        )

        config_data = results[
            encoding
        ][
            config_name
        ]

        if config_data.get("encoding_error"):
            continue

        heatmap = np.array([
            [
                config_data[
                    "profiles"
                ][
                    profile_name
                ][
                    error_rate
                ][
                    "recovery_rate"
                ]
                for error_rate in ERROR_RATES
            ]
            for profile_name in ERROR_PROFILES
        ])

        plt.figure(figsize=(16, 6))

        image = plt.imshow(
            heatmap,
            aspect="auto",
            interpolation="nearest",
            vmin=0,
            vmax=1,
        )

        plt.xticks(
            range(len(ERROR_RATES)),
            [f"{rate:g}" for rate in ERROR_RATES],
            rotation=60,
            ha="right",
        )

        plt.yticks(
            range(len(ERROR_PROFILES)),
            list(ERROR_PROFILES.keys())
        )

        plt.xlabel(
            "Total sequencing error rate"
        )

        plt.ylabel(
            "Error composition"
        )

        plt.title(
            f"E2: {encoding} — recovery heatmap\n"
            f"Reference configuration: {config_name}, mean={FIXED_MEAN}"
        )

        colorbar = plt.colorbar(image)
        colorbar.set_label("Recovery rate")

        plt.tight_layout()

        png_file = os.path.join(
            plot_dir,
            f"E2_{encoding}_heatmap.png"
        )

        pdf_file = os.path.join(
            plot_dir,
            f"E2_{encoding}_heatmap.pdf"
        )

        plt.savefig(
            png_file,
            dpi=300,
            bbox_inches="tight"
        )

        plt.savefig(
            pdf_file,
            bbox_inches="tight"
        )

        plt.close()

        print(f"Saved: {png_file}")
        print(f"Saved: {pdf_file}")


# ================================================================
# PLOT 3 — SCATTER PLOT
# ================================================================
#
# Numeric relationship:
#   total error rate vs recovery rate
#
# Each point is a profile/error-rate observation from the reference
# configuration.
# ================================================================


def plot_reference_scatter(results):

    plot_dir = os.path.join(
        OUTPUT_DIR,
        "plots"
    )

    os.makedirs(plot_dir, exist_ok=True)

    plt.figure(figsize=(12, 8))

    for encoding in ENCODINGS:

        config_name = first_configuration_name(
            results,
            encoding
        )

        config_data = results[
            encoding
        ][
            config_name
        ]

        if config_data.get("encoding_error"):
            continue

        for profile_name in ERROR_PROFILES:

            profile_data = config_data[
                "profiles"
            ][
                profile_name
            ]

            x = []
            y = []

            for error_rate, point in profile_data.items():

                x.append(error_rate)
                y.append(point["recovery_rate"])

            # Use encoding as the marker label and profile as point
            # transparency/shape information without forcing 49 legend
            # entries.
            plt.scatter(
                x,
                y,
                s=35,
                alpha=0.55,
                label=encoding if profile_name == "Balanced" else None,
            )

    plt.xlabel(
        "Total sequencing error rate"
    )

    plt.ylabel(
        "Recovery rate"
    )

    plt.title(
        "E2: Recovery scatter plot across encodings and error compositions\n"
        f"Reference configurations, mean={FIXED_MEAN}"
    )

    plt.xscale("log")
    plt.ylim(0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend(title="Encoding")
    plt.tight_layout()

    png_file = os.path.join(
        plot_dir,
        "E2_reference_scatter.png"
    )

    pdf_file = os.path.join(
        plot_dir,
        "E2_reference_scatter.pdf"
    )

    plt.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")


# ================================================================
# PLOT 4 — BAR PLOT AT 1% TOTAL ERROR
# ================================================================


def plot_one_percent_bar(results):

    comparison_error_rate = 0.01

    plot_dir = os.path.join(
        OUTPUT_DIR,
        "plots"
    )

    os.makedirs(plot_dir, exist_ok=True)

    profile_names = list(
        ERROR_PROFILES.keys()
    )

    encoding_names = list(
        ENCODINGS
    )

    n_profiles = len(profile_names)
    n_encodings = len(encoding_names)

    x = np.arange(n_profiles)

    width = 0.8 / n_encodings

    plt.figure(figsize=(15, 7))

    for encoding_index, encoding in enumerate(encoding_names):

        config_name = first_configuration_name(
            results,
            encoding
        )

        config_data = results[
            encoding
        ][
            config_name
        ]

        values = []

        for profile_name in profile_names:

            point = config_data[
                "profiles"
            ][
                profile_name
            ].get(comparison_error_rate)

            values.append(
                point["recovery_rate"]
                if point is not None
                else 0.0
            )

        offset = (
            encoding_index
            - (n_encodings - 1) / 2
        ) * width

        plt.bar(
            x + offset,
            values,
            width=width,
            label=encoding,
        )

    plt.xticks(
        x,
        profile_names,
        rotation=25,
        ha="right"
    )

    plt.xlabel(
        "Error composition"
    )

    plt.ylabel(
        "Recovery rate"
    )

    plt.title(
        "E2: Recovery at 1% total sequencing error\n"
        f"mean={FIXED_MEAN}"
    )

    plt.ylim(0, 1.05)
    plt.grid(axis="y", alpha=0.3)
    plt.legend(title="Encoding")
    plt.tight_layout()

    png_file = os.path.join(
        plot_dir,
        "E2_recovery_at_1pct_bar.png"
    )

    pdf_file = os.path.join(
        plot_dir,
        "E2_recovery_at_1pct_bar.pdf"
    )

    plt.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")


# ================================================================
# COST VS ERROR TOLERANCE
# ================================================================
#
# Similar to E1, but now tolerance is tracked separately for each
# error composition profile.
# ================================================================


def calculate_cost_vs_tolerance(results, target_recovery):

    output = {}

    for encoding, encoding_data in results.items():

        output[encoding] = {}

        for config_name, config_data in encoding_data.items():

            if config_data.get("encoding_error"):
                continue

            dna_cost = config_data[
                "dna_cost_nt_per_bit"
            ]

            output[encoding][config_name] = {}

            for profile_name, profile_data in config_data[
                "profiles"
            ].items():

                tolerated_error_rate = None
                success_at_tolerance = None

                for error_rate in sorted(profile_data.keys()):

                    point = profile_data[
                        error_rate
                    ]

                    success_rate = point[
                        "recovery_rate"
                    ]

                    if success_rate >= target_recovery:

                        tolerated_error_rate = error_rate
                        success_at_tolerance = success_rate

                    else:
                        # Error rates are increasing, matching E1's
                        # tolerance calculation pattern.
                        break

                output[encoding][config_name][profile_name] = {
                    "dna_cost": dna_cost,
                    "error_tolerance": tolerated_error_rate,
                    "success_rate": success_at_tolerance,
                }

    return output


# ================================================================
# PLOT COST VS ERROR TOLERANCE
# ================================================================


def plot_cost_vs_tolerance(
    cost_vs_tolerance,
    target_recovery,
):

    plot_dir = os.path.join(
        OUTPUT_DIR,
        "plots"
    )

    os.makedirs(plot_dir, exist_ok=True)

    plt.figure(figsize=(15, 9))

    for encoding, config_data in cost_vs_tolerance.items():

        for config_name, profile_data in config_data.items():

            for profile_name, result in profile_data.items():

                tolerance = result[
                    "error_tolerance"
                ]

                dna_cost = result[
                    "dna_cost"
                ]

                if tolerance is None:
                    continue

                # Plot all configurations/profile combinations.
                # The label is kept compact; raw CSV contains the
                # complete configuration information.
                plt.scatter(
                    dna_cost,
                    tolerance,
                    s=55,
                    alpha=0.7,
                    label=(
                        encoding
                        if profile_name == "Balanced"
                        else None
                    ),
                )

    plt.xlabel(
        "DNA cost (nt/bit)"
    )

    plt.ylabel(
        "Maximum tolerated total error rate"
    )

    plt.title(
        f"E2: DNA cost vs error tolerance\n"
        f"Target recovery = {target_recovery:.0%}"
    )

    plt.grid(True, alpha=0.3)
    plt.legend(title="Encoding")
    plt.tight_layout()

    recovery_percent = int(
        target_recovery * 100
    )

    png_file = os.path.join(
        plot_dir,
        f"E2_cost_vs_tolerance_{recovery_percent}pct.png"
    )

    pdf_file = os.path.join(
        plot_dir,
        f"E2_cost_vs_tolerance_{recovery_percent}pct.pdf"
    )

    plt.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")


# ================================================================
# SAVE COST VS TOLERANCE CSV
# ================================================================


def save_cost_vs_tolerance_csv(cost_vs_tolerance):

    output_file = os.path.join(
        OUTPUT_DIR,
        f"E2_cost_vs_tolerance_{TIMESTAMP}.csv"
    )

    with open(
        output_file,
        "w",
        newline="",
    ) as f:

        writer = csv.writer(f)

        writer.writerow([
            "target_recovery",
            "encoding",
            "configuration",
            "error_profile",
            "dna_cost_nt_per_bit",
            "error_tolerance",
            "success_rate_at_tolerance",
        ])

        for target_recovery, target_data in cost_vs_tolerance.items():

            for encoding, encoding_data in target_data.items():

                for config_name, profile_data in encoding_data.items():

                    for profile_name, result in profile_data.items():

                        writer.writerow([
                            target_recovery,
                            encoding,
                            config_name,
                            profile_name,
                            result["dna_cost"],
                            result["error_tolerance"],
                            result["success_rate"],
                        ])

    print(
        f"Saved cost/tolerance CSV: {output_file}"
    )


# ================================================================
# MAIN
# ================================================================


if __name__ == "__main__":

    print()
    print("=" * 90)
    print("STARTING E2")
    print("=" * 90)
    print()

    print(
        f"Input file: {INPUT_FILENAME}"
    )

    if not os.path.isfile(
        os.path.join(
            "./tests/testfiles",
            INPUT_FILENAME
        )
    ):
        raise FileNotFoundError(
            f"Could not find input file: "
            f"./tests/testfiles/{INPUT_FILENAME}"
        )

    print()

    results = simulate_e2()

    print_summary(results)

    # ------------------------------------------------------------
    # Save raw simulation data.
    # ------------------------------------------------------------

    save_json(results)

    save_summary_csv(results)

    # ------------------------------------------------------------
    # Main E2 plots.
    # ------------------------------------------------------------

    plot_reference_recovery_curves(results)

    plot_reference_heatmaps(results)

    plot_reference_scatter(results)

    plot_one_percent_bar(results)

    # ------------------------------------------------------------
    # E1-style cost vs error-tolerance analysis.
    # ------------------------------------------------------------

    cost_vs_tolerance = {
        1.0: calculate_cost_vs_tolerance(
            results,
            1.0
        ),
        0.90: calculate_cost_vs_tolerance(
            results,
            0.90
        ),
    }

    plot_cost_vs_tolerance(
        cost_vs_tolerance[1.0],
        1.0,
    )

    plot_cost_vs_tolerance(
        cost_vs_tolerance[0.90],
        0.90,
    )

    save_cost_vs_tolerance_csv(
        cost_vs_tolerance
    )

    print()
    print("=" * 90)
    print("E2 COMPLETE")
    print("=" * 90)
    print()
    print(f"Results saved in: {OUTPUT_DIR}")
