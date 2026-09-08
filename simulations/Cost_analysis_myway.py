
import json
import os
import csv
import copy
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.encode import Encode
from dnabyte.params import Params
from dnabyte.sequence import SimulateSequencing
from dnabyte.store import SimulateStorage
from dnabyte.synthesize import SimulateSynthesis
from simulations.simulation import Simulation
from dnabyte.binarize import Binarize
from dnabyte.data_classes import Data
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus

def encode_and_count_single_encoding(encoding, defultparams):
    # Sanity check for encoding method
    sim = Simulation([defultparams[encoding]])
    results = sim.run()
    
    # Check if any simulation succeeded
    encsuccess = any(
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
    encodedbp = sum(len(i) for i in data_enc.data)*defultparams[encoding].mean  # Sum lengths of all encoded sequences
    encodedbpnomean = sum(len(i) for i in data_enc.data)  # Sum lengths of all encoded sequences without mean

    print(f"Encoded bp count for {encoding}: {encodedbp}")
    print(f"Encoding success status for {encoding}: {encsuccess}")

    if encsuccess:
        return encodedbp, encodedbpnomean   
    else:
        return None  # Return None if the encoding failed

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
        encodedbp[encoding] = sum(len(i) for i in data_enc.data)*defultparams[encoding].mean  # Sum lengths of all encoded sequences

    print(f"Encoded bp counts: {encodedbp}")
    print(f"Encoding success status: {encsuccess}")

    if all(value is True for value in encsuccess.values()):
        return encodedbp
    else:
        failed_encodings = [encoding for encoding, success in encsuccess.items() if not success]
        raise RuntimeError(f"Baseline simulation failed for the following encodings: {', '.join(failed_encodings)}")

def find_smallest_target(d, tolerance=0.1, max_target=1000):
    start = max(d.values())  # multiplier must be >= 1

    for target in range(start, max_target + 1):
        multipliers = {}
        scaled_values = {}

        for name, value in d.items():
            multiplier = round(target / value)  # integer
            scaled = value * multiplier

            multipliers[name] = multiplier
            scaled_values[name] = scaled

        # Check whether all scaled values are sufficiently close
        max_error = max(
            abs(scaled - target) / target
            for scaled in scaled_values.values()
        )

        if max_error <= tolerance:
            return {
                "target": target,
                "multipliers": multipliers,
                "scaled_values": scaled_values,
                "max_error": max_error,
            }

    return None

    return None
#Function with the found amount of bp that increses the enecoding parameters of all other encodings to have the same amount of bp and then run all encodings with error correction and find the highest strand/bp count;
def encode_and_count_with_error_correction(encodings, defultparams, noECencodinglengths):
    final_encodedbp = {}
    final_encodedbpnomean = {}
    maxbp = 10000
    target_bp = max(noECencodinglengths.values())
    max_encoding = max(noECencodinglengths, key=noECencodinglengths.get)
    final_encodedbp[max_encoding] = target_bp
    final_encodedbpnomean[max_encoding] = target_bp
    for encoding in encodings:
        if encoding == 'goldman' or encoding == 'yinyang':
            if encoding == 'yinyang':
                defultparams[encoding].mean = 1
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)

            previous_success = success
            previous_encodedbpnomean = encodedbpnomean
            
            if success is None:
                raise RuntimeError(f"Encoding failed for {encoding}")
            # while True:
            #     if success < target_bp+300:
            #         previous_success = success
            #         previous_encodedbpnomean = encodedbpnomean
            #         defultparams[encoding].mean += 1
            #     elif success > target_bp+300:
            #         if abs(previous_success - target_bp) < abs(success - target_bp):
            #             defultparams[encoding].mean -= 1
            #             success = previous_success
            #         break
            #     success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            #     if success is None:
            #         raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean
            
        elif encoding == 'church' or encoding == 'wukong':
            defultparams[encoding].add_redundancy = True
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding add_redundancy failed for {encoding}")
            # if success > target_bp+300:
            #     defultparams[encoding].add_redundancy = False
            #     success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            #     if success is None:
            #         raise RuntimeError(f"Encoding add_redundancy failed for {encoding}")
            while True:
                
                defultparams[encoding].rs_num += 1
                
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    defultparams[encoding].rs_num -= 1
                    break
                if maxbp < success:
                    break   
                    # raise RuntimeError(f"Encoding add_redundancy failed for {encoding}")

            # while True:
            #     if success < target_bp+300:
            #         previous_success = success
            #         defultparams[encoding].mean += 1
            #     elif success > target_bp+300:
            #         if abs(previous_success - target_bp) < abs(success - target_bp):
            #             defultparams[encoding].mean -= 1
            #             success = previous_success
            #         break
            #     success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            #     if success is None:
            #         raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean

        elif encoding == 'gcplus':
            defultparams[encoding].gcplus_c1 = 1
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding gcplus_c1 failed for {encoding}")
            while True:
                # if success < target_bp+300:
                #     previous_success = success
                #     defultparams[encoding].gcplus_c1 += 1
                # elif success > target_bp+300:
                #     if abs(previous_success - target_bp) < abs(success - target_bp):
                defultparams[encoding].gcplus_c1 += 1
                    #     success = previous_success
                    # break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    defultparams[encoding].gcplus_c1 -= 1
                    break
                if maxbp < success:
                    break  
                    # raise RuntimeError(f"Encoding gcplus_c1 failed for {encoding}")
            # while True:
            #     if success < target_bp+300:
            #         previous_success = success
            #         defultparams[encoding].mean += 1
            #     elif success > target_bp+300:
            #         if abs(previous_success - target_bp) < abs(success - target_bp):
            #             defultparams[encoding].mean -= 1
            #             success = previous_success
            #         break
            #     success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            #     if success is None:
            #         raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean

        elif encoding == 'max_density' or encoding == 'no_homopolymer':
            defultparams[encoding].clustering_method = None
            defultparams[encoding].recovery_method = None
            defultparams[encoding].outer_error_correction = 'reedsolomon'
            defultparams[encoding].reed_solo_percentage = 0.8
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding max_density failed for {encoding}")
            while True:
                # if success < target_bp+300:
                #     previous_success = success
                defultparams[encoding].reed_solo_percentage = defultparams[encoding].reed_solo_percentage - 0.01
                # elif success > target_bp+300:
                #     if abs(previous_success - target_bp) < abs(success - target_bp):
                #         defultparams[encoding].reed_solo_percentage = min(defultparams[encoding].reed_solo_percentage + 0.01, 1.0)
                #         success = previous_success  
                #     break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    defultparams[encoding].reed_solo_percentage = defultparams[encoding].reed_solo_percentage + 0.01
                    break
                if maxbp < success:
                    break  
                    # raise RuntimeError(f"Encoding failed for {encoding}")
            # while True:
            #     if success < target_bp+300:
            #         previous_success = success
            #         defultparams[encoding].mean += 1
            #     elif success > target_bp+300:
            #         if abs(previous_success - target_bp) < abs(success - target_bp):
            #             defultparams[encoding].mean -= 1
            #             success = previous_success
            #         break
            #     success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            #     if success is None:
            #         raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean

        elif encoding == 'hedeges':
            defultparams[encoding].hedges_coderate = 3
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding hedges failed for {encoding}")
            while True:
                # previous_success = success
                # if success < target_bp+300:
                defultparams[encoding].hedges_coderate += 1
                # elif success > target_bp+300:
                #     if abs(previous_success - target_bp) < abs(success - target_bp):
                #         defultparams[encoding].hedges_coderate -= 1
                #         success = previous_success
                #     break

                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    defultparams[encoding].hedges_coderate -= 1
                    break
                if maxbp < success:
                    break  
                    # raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean

    print(f"Final encoded bp counts with error correction: {final_encodedbp}")
    print(f"Final encoded bp counts without mean: {final_encodedbpnomean}")

    return defultparams
# main function run for the found highest strand/bp count for diffrent error rates and find the highest strand/bp count for each encoding with error correction and without error correction;
def simulate_cost_analysis(encodings, defultparams, error_rates):
    error_simulation_results = {}
    for error_rate in error_rates:
        print(f"Simulating for error rate: {error_rate}")
        for encoding in encodings:
            defultparams[encoding].iid_error_rate = error_rate
            defultparams[encoding].iid_substitution_rate = error_rate
            defultparams[encoding].iid_insertion_rate = 0.0
            defultparams[encoding].iid_deletion_rate = 0.0

            if defultparams[encoding].mean > 1:
                defultparams[encoding].similarity_threshold = 0.75

            binarizer = Binarize(defultparams[encoding])

            if hasattr(defultparams[encoding], 'file_paths') and defultparams[encoding].file_paths:
                file_paths = ['./tests/testfiles/' + fp for fp in defultparams[encoding].file_paths]
            elif hasattr(defultparams[encoding], 'filename'):
                file_paths = ['./tests/testfiles/' + defultparams[encoding].filename]
            else:
                raise ValueError("params must have either 'file_paths' or 'filename' attribute")

            if encoding not in error_simulation_results:
                error_simulation_results[encoding] = {}
            if error_rate not in error_simulation_results[encoding]:
                error_simulation_results[encoding][error_rate] = []

            data_obj = Data(file_paths=file_paths)
            binary_code = binarizer.binarize(data_obj)
            coder = Encode(defultparams[encoding])
            data_enc, info = coder.encode(binary_code)
            information_bits = len(binary_code.data)

            # --------------------------------------------------
            # CALCULATE DNA STORAGE COST
            # --------------------------------------------------

            dna_length = sum(
                len(str(sequence))
                for sequence in data_enc
            )

            dna_cost_nt_per_bit = (
                dna_length / information_bits
                if information_bits > 0
                else 0
            )

            # Expected sequencing errors for this encoded dataset
            expected_errors = dna_length * error_rate

            for i in range(50):
                success = False
                res = None

                try:
                    syn = SimulateSynthesis(defultparams[encoding])
                    data_syn, info = syn.simulate(data_enc)

                    sto = SimulateStorage(defultparams[encoding])
                    data_sto, info = sto.simulate(data_syn)

                    seq = SimulateSequencing(defultparams[encoding])
                    data_seq, info = seq.simulate(data_sto)

                    if defultparams[encoding].clustering_method is not None and defultparams[encoding].recovery_method is not None:
                        cluster_obj = Cluster(defultparams[encoding])
                        data_cluster, info = cluster_obj.cluster(data_seq)
                        consensus_obj = Consensus(defultparams[encoding])
                        data_cor, info = consensus_obj.call(data_cluster)
                    else:
                        # Fall back to encoding-specific process
                        data_cor, info = coder.process(data_seq)

                    data_dec, valid, info = coder.decode(data_cor)
                    if not valid:
                        raise ValueError("decoded data invalid")

                    comparison, res = data_dec.compare(data_dec, binary_code)
                    success = comparison == 'SUCCESS'
                    # print(binary_code.data)
                    # print(data_dec.data)

                except Exception as e:
                    import traceback
                    print(f"\n❌ Run {i+1} failed for {encoding} at error rate {error_rate}")
                    print(f"Error: {type(e).__name__}: {e}")
                    print(f"Traceback: {traceback.format_exc()}\n")
                    success = False
                    res = None

                error_simulation_results[encoding][error_rate].append({
                    "run": i + 1,
                    "success": success,
                    "results": res,

                    # Normalization information
                    "information_bits": information_bits,
                    "dna_length": dna_length,
                    "dna_cost_nt_per_bit": dna_cost_nt_per_bit,

                    # Expected number of nucleotide errors
                    "expected_errors": expected_errors,
                })
                


            # for i in range(50):  # Run each simulation 20 times for averaging
            #     print(f"Run {i+1} for encoding {encoding} at error rate {error_rate}")
            #     sim = Simulation([defultparams[encoding]])
            #     results = sim.run()
                
            #     # Check if any simulation succeeded
            #     encsuccess = any(
            #         sim_result.get('status') == 'SUCCESS' 
            #         for sim_result in results.values()
            #     )
                
            #     if encoding not in error_simulation_results:
            #         error_simulation_results[encoding] = {}
                
            #     if error_rate not in error_simulation_results[encoding]:
            #         error_simulation_results[encoding][error_rate] = []
                
            #     error_simulation_results[encoding][error_rate].append({
            #         "run": i + 1,
            #         "success": encsuccess,
            #         "results": results
            #     })
            
            # # Check if any simulation succeeded
            # encsuccess = any(
            #     sim_result.get('status') == 'SUCCESS' 
            #     for sim_result in results.values()
            # )

    return error_simulation_results

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
            'kmer_size_cluster': 180,
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
            'gcplus_c1': 2,
            'barcode_length': 0,
            'left_primer': '',
            'right_primer': '',
            'kmer_size_debruijn': 90,
            # 'kmer_size_cluster': 180,
            'kmer_size_cluster': 20,
            'kmer_threshold': 0.9,
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
def get_parameter_configurations(encoding):

    if encoding == 'church':
        return [
            {
                'add_redundancy': add_redundancy,
                'rs_num': rs_num
            }
            for add_redundancy in [False, True]
            for rs_num in range(0, 6)
        ]

    elif encoding == 'wukong':
        return [
            {
                'add_redundancy': add_redundancy,
                'rs_num': rs_num
            }
            for add_redundancy in [False, True]
            for rs_num in range(0, 6)
        ]

    elif encoding == 'goldman':
        return [
            {
                'mean': mean
            }
            for mean in [1]
        ]

    elif encoding == 'gcplus':
        return [
            {
                'gcplus_c1': c1
            }
            for c1 in [1, 2, 3, 4, 5]
        ]

    elif encoding == 'hedges':
        return [
            {
                'hedges_coderate': rate
            }
            for rate in [1, 2, 3, 4, 5, 6]
        ]

    elif encoding in ['max_density', 'no_homopolymer']:
        return [
            {
                'reed_solo_percentage': percentage
            }
            for percentage in [
                0.50, 0.60, 0.70, 0.80, 0.90
            ]
        ]

    else:
        raise ValueError(
            f"No parameter sweep defined for {encoding}"
        )



def simulate_cost_analysis(encodings, base_params, error_rates):
    """
    Test every parameter configuration for every encoding at every
    sequencing error rate.

    Result structure:

        results[encoding][configuration][error_rate] = [
            run1,
            run2,
            ...
        ]

    The DNA cost is measured from the encoded DNA BEFORE sequencing
    errors are introduced.
    """

    error_simulation_results = {}

    for encoding in encodings:

        print("\n" + "=" * 80)
        print(f"ENCODING: {encoding}")
        print("=" * 80)

        error_simulation_results[encoding] = {}

        # ------------------------------------------------------------
        # Get all parameter combinations for this encoding
        # ------------------------------------------------------------

        parameter_configs = get_parameter_configurations(encoding)

        for config in parameter_configs:

            # --------------------------------------------------------
            # Create a completely independent copy of the parameters
            # --------------------------------------------------------

            params = copy.deepcopy(base_params[encoding])

            # Apply the parameters for this configuration
            for parameter, value in config.items():
                setattr(params, parameter, value)

            # --------------------------------------------------------
            # Give this configuration a readable name
            # --------------------------------------------------------

            config_name = "_".join(
                f"{key}={value}"
                for key, value in config.items()
            )

            print("\n" + "-" * 80)
            print(f"Configuration: {config_name}")
            print("-" * 80)

            error_simulation_results[encoding][config_name] = {}

            # --------------------------------------------------------
            # IMPORTANT:
            #
            # For max_density and no_homopolymer, changing
            # reed_solo_percentage does nothing unless Reed-Solomon
            # outer error correction is actually enabled.
            # --------------------------------------------------------

            if encoding in ['max_density', 'no_homopolymer']:

                params.outer_error_correction = 'reedsolomon'

                # These encodings do not use the normal clustering /
                # recovery pipeline in your experiment.
                params.clustering_method = None
                params.recovery_method = None

            # --------------------------------------------------------
            # Test every sequencing error rate
            # --------------------------------------------------------

            for error_rate in error_rates:

                print(
                    f"  Error rate: {error_rate}"
                )

                # ----------------------------------------------------
                # Set IID sequencing error parameters
                # ----------------------------------------------------

                params.iid_error_rate = error_rate
                params.iid_substitution_rate = error_rate
                params.iid_insertion_rate = 0.0
                params.iid_deletion_rate = 0.0

                # ----------------------------------------------------
                # Adjust similarity threshold for encodings that
                # use mean > 1
                #
                # getattr() is used because not every encoding
                # necessarily has a 'mean' parameter.
                # ----------------------------------------------------

                if getattr(params, 'mean', 1) > 1:
                    params.similarity_threshold = 0.75

                # ----------------------------------------------------
                # BINARIZE INPUT FILE
                # ----------------------------------------------------

                binarizer = Binarize(params)

                if (
                    hasattr(params, 'file_paths')
                    and params.file_paths
                ):
                    file_paths = [
                        './tests/testfiles/' + fp
                        for fp in params.file_paths
                    ]

                elif hasattr(params, 'filename'):
                    file_paths = [
                        './tests/testfiles/' + params.filename
                    ]

                else:
                    raise ValueError(
                        "params must have either "
                        "'file_paths' or 'filename' attribute"
                    )

                data_obj = Data(file_paths=file_paths)

                binary_code = binarizer.binarize(data_obj)

                # ----------------------------------------------------
                # ENCODE
                # ----------------------------------------------------

                coder = Encode(params)

                data_enc, info = coder.encode(binary_code)

                # Number of original information bits
                information_bits = len(binary_code.data)

                # ----------------------------------------------------
                # DNA LENGTH
                #
                # IMPORTANT:
                #
                # This must be measured BEFORE sequencing errors.
                # Otherwise insertions/deletions would incorrectly
                # change the storage-cost measurement.
                # ----------------------------------------------------

                if hasattr(data_enc, 'data'):
                    dna_sequences = data_enc.data
                else:
                    dna_sequences = data_enc

                dna_length = sum(
                    len(str(sequence))
                    for sequence in dna_sequences
                )

                # ----------------------------------------------------
                # DNA STORAGE COST
                #
                # nt per original information bit
                # ----------------------------------------------------

                dna_cost_nt_per_bit = (
                    dna_length / information_bits
                    if information_bits > 0
                    else 0
                )

                # ----------------------------------------------------
                # Expected number of sequencing errors across all
                # synthesized nucleotides
                # ----------------------------------------------------

                expected_errors = dna_length * error_rate

                results_for_runs = []

                # ----------------------------------------------------
                # Repeat the experiment multiple times because the
                # IID sequencing errors are random.
                # ----------------------------------------------------

                for i in range(50):

                    success = False
                    res = None

                    try:

                        # ==================================================
                        # SYNTHESIS
                        # ==================================================

                        syn = SimulateSynthesis(params)

                        data_syn, info = syn.simulate(data_enc)

                        # ==================================================
                        # STORAGE
                        # ==================================================

                        sto = SimulateStorage(params)

                        data_sto, info = sto.simulate(data_syn)

                        # ==================================================
                        # SEQUENCING
                        # ==================================================

                        seq = SimulateSequencing(params)

                        data_seq, info = seq.simulate(data_sto)

                        # ==================================================
                        # RECOVERY / CLUSTERING
                        # ==================================================

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

                            data_cor, info = coder.process(data_seq)

                        # ==================================================
                        # DECODE
                        # ==================================================

                        data_dec, valid, info = coder.decode(data_cor)

                        if not valid:
                            raise ValueError(
                                "decoded data invalid"
                            )

                        # ==================================================
                        # COMPARE DECODED DATA WITH ORIGINAL DATA
                        # ==================================================

                        comparison, res = data_dec.compare(
                            data_dec,
                            binary_code
                        )

                        success = comparison == 'SUCCESS'

                    except Exception as e:

                        print(
                            f"    Run {i + 1} failed "
                            f"for {encoding}, "
                            f"{config_name}, "
                            f"error rate {error_rate}: "
                            f"{type(e).__name__}: {e}"
                        )

                        success = False
                        res = None

                    # ------------------------------------------------
                    # Store result of this run
                    # ------------------------------------------------

                    results_for_runs.append({
                        "run": i + 1,
                        "success": success,
                        "results": res,

                        "information_bits": information_bits,

                        "dna_length": dna_length,

                        "dna_cost_nt_per_bit":
                            dna_cost_nt_per_bit,

                        "expected_errors":
                            expected_errors,

                        "error_rate":
                            error_rate,

                        "encoding":
                            encoding,

                        "configuration":
                            config_name,
                    })

                # ----------------------------------------------------
                # Store all 50 runs for this error rate
                # ----------------------------------------------------

                error_simulation_results[
                    encoding
                ][
                    config_name
                ][
                    error_rate
                ] = results_for_runs

    return error_simulation_results
#, 'hedges', 'wukong', 'no_homopolymer', 'church', 'max_density', 'goldman'

if __name__ == '__main__':

    # ================================================================
    # ENCODINGS TO TEST
    # ================================================================

    encodings = [
        'goldman',
        'church',
        'no_homopolymer'
    ]

    # ================================================================
    # COMMON PARAMETERS
    # ================================================================

    base_params = {
        'filename': 'testfilesimsall.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'iid',
        'recovery_method': 'pass_through',
        'min_coverage': 1,
        'clustering_method': 'pass_through',
        'storage_conditions': None,
        'kmer_seed': 42,
        'similarity_threshold': 1,
        'gap_penalty': 1.0,
    }

    # ================================================================
    # DEFAULT ERROR / REPLICATION PARAMETERS
    # ================================================================

    no_error_params = {
        'mean': 1,
        'std_dev': 0,

        'iid_error_rate': 0.0,
        'iid_substitution_rate': 0.0,
        'iid_insertion_rate': 0.0,
        'iid_deletion_rate': 0.0,
    }

    # ================================================================
    # SEQUENCING ERROR RATES
    # ================================================================

    error_rates = [
        0.0,
        0.0001,
        0.0005,
        0.001,
        0.005,
        0.01
    ]

    # ================================================================
    # BUILD PARAMETERS FOR EACH ENCODING
    # ================================================================

    defoultparams = {}

    for encoding in encodings:

        parameters = base_params.copy()

        parameters.update(
            default_values_of_encoding(encoding)
        )

        parameters.update(
            no_error_params
        )

        # ------------------------------------------------------------
        # max_density and no_homopolymer use their own outer
        # error-correction / decoding pipeline.
        #
        # The actual Reed-Solomon settings are applied inside
        # simulate_cost_analysis() for every sweep configuration.
        # ------------------------------------------------------------

        if encoding in [
            'max_density',
            'no_homopolymer'
        ]:

            parameters['clustering_method'] = None
            parameters['recovery_method'] = None

        print(
            f"\nParameters for {encoding}:"
        )

        print(parameters)

        defoultparams[encoding] = Params(
            **parameters
        )

    # ================================================================
    # RUN THE COMPLETE PARAMETER SWEEP
    # ================================================================
    #
    # IMPORTANT:
    #
    # Do NOT call encode_and_count_with_error_correction() here.
    #
    # That function mutates defoultparams and would mean that the
    # parameter sweep no longer starts from the original baseline
    # configuration.
    #
    # simulate_cost_analysis() makes its own deepcopy for every
    # configuration.
    # ================================================================

    sim_resoults = simulate_cost_analysis(
        encodings,
        defoultparams,
        error_rates
    )

    print("\n")
    print("=" * 80)
    print("PARAMETER SWEEP FINISHED")
    print("=" * 80)
    # ================================================================
    # SAVE RAW SIMULATION RESULTS
    # ================================================================
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "simulation_results_costs")
    os.makedirs(
        output_dir,
        exist_ok=True
    )

    raw_results_file = os.path.join(
        output_dir,
        "simulation_results.json"
    )

    with open(raw_results_file, "w") as f:
        json.dump(
            sim_resoults,
            f,
            indent=4
        )

    print(f"Saved raw results: {raw_results_file}")

    # ================================================================
    # PRINT RESULTS
    # ================================================================

    for encoding, config_data in sim_resoults.items():

        print("\n")
        print("#" * 80)
        print(f"ENCODING: {encoding}")
        print("#" * 80)

        for config_name, error_data in config_data.items():

            print(
                f"\nConfiguration: {config_name}"
            )

            for error_rate, runs in error_data.items():

                successful_runs = sum(
                    run['success']
                    for run in runs
                )

                total_runs = len(runs)

                success_rate = (
                    successful_runs / total_runs
                    if total_runs > 0
                    else 0
                )

                dna_cost = runs[0][
                    'dna_cost_nt_per_bit'
                ]

                print(
                    f"  Error rate: {error_rate:<8} "
                    f"Success: "
                    f"{successful_runs}/{total_runs} "
                    f"({success_rate:.2%}) "
                    f"DNA cost: "
                    f"{dna_cost:.4f} nt/bit"
                )
    import matplotlib.pyplot as plt

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(script_dir, "simulation_results_costs")
    os.makedirs(output_dir, exist_ok=True)

    # ================================================================
    # SUCCESS RATE VS ERROR RATE
    # ================================================================

    plt.figure(figsize=(14, 9))

    for encoding, config_data in sim_resoults.items():

        for config_name, error_data in config_data.items():

            sorted_error_rates = sorted(error_data.keys())

            success_rates = []

            for error_rate in sorted_error_rates:

                runs = error_data[error_rate]

                successful_runs = sum(
                    run['success']
                    for run in runs
                )

                total_runs = len(runs)

                success_rate = (
                    successful_runs / total_runs
                    if total_runs > 0
                    else 0
                )

                success_rates.append(success_rate)

            plt.plot(
                sorted_error_rates,
                success_rates,
                'o-',
                label=f"{encoding}: {config_name}"
            )


    plt.xlabel(
        "IID substitution error rate"
    )

    plt.ylabel(
        "Successful recovery rate"
    )

    plt.title(
        "Recovery success rate vs sequencing error rate"
    )

    plt.ylim(
        0,
        1.05
    )

    plt.grid(
        True,
        alpha=0.3
    )

    plt.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left'
    )

    plt.tight_layout()


    # ================================================================
    # SAVE FIGURE
    # ================================================================

    # Make sure the output directory exists
    os.makedirs(
        output_dir,
        exist_ok=True
    )

    # PNG path
    success_png = os.path.join(
        output_dir,
        "success_rate_vs_error_rate.png"
    )

    # PDF path
    success_pdf = os.path.join(
        output_dir,
        "success_rate_vs_error_rate.pdf"
    )


    # SAVE BEFORE plt.show()
    plt.savefig(
        success_png,
        dpi=300,
        bbox_inches='tight'
    )

    plt.savefig(
        success_pdf,
        bbox_inches='tight'
    )


    # Confirm that the files actually exist
    print()
    print("================================================")
    print("SUCCESS RATE PLOT SAVED")
    print("================================================")
    print(f"PNG: {success_png}")
    print(f"PDF: {success_pdf}")
    print(f"PNG exists: {os.path.exists(success_png)}")
    print(f"PDF exists: {os.path.exists(success_pdf)}")
    print("================================================")


    # Show the figure AFTER saving
    plt.show()

    # Close after showing
    plt.close()

    # ================================================================
    # DNA COST VS ERROR TOLERANCE
    # ================================================================

    target_recovery = 1

    cost_vs_tolerance = {}


    for encoding, config_data in sim_resoults.items():

        cost_vs_tolerance[encoding] = {}

        for config_name, error_data in config_data.items():

            # --------------------------------------------------------
            # DNA cost is independent of sequencing error rate, so
            # take it from the first available result.
            # --------------------------------------------------------

            first_error_rate = next(
                iter(error_data)
            )

            first_run = error_data[first_error_rate][0]

            dna_cost = first_run[
                'dna_cost_nt_per_bit'
            ]

            # --------------------------------------------------------
            # Find the highest tested error rate for which the
            # configuration still achieves the target recovery rate.
            # --------------------------------------------------------

            sorted_error_rates = sorted(
                error_data.keys()
            )

            tolerated_error_rate = None
            tolerated_success_rate = None

            for error_rate in sorted_error_rates:

                runs = error_data[error_rate]

                successful_runs = sum(
                    run['success']
                    for run in runs
                )

                total_runs = len(runs)

                success_rate = (
                    successful_runs / total_runs
                    if total_runs > 0
                    else 0
                )

                if success_rate >= target_recovery:

                    tolerated_error_rate = error_rate
                    tolerated_success_rate = success_rate

                else:

                    # Once the success rate drops below the target,
                    # stop because the error rates are increasing.
                    break

            cost_vs_tolerance[encoding][config_name] = {
                'dna_cost': dna_cost,
                'error_tolerance': tolerated_error_rate,
                'success_rate': tolerated_success_rate
            }

    # ================================================================
    # SAVE COST VS ERROR TOLERANCE DATA
    # ================================================================

    summary_file = os.path.join(
        output_dir,
        "cost_vs_error_tolerance.csv"
    )

    with open(
        summary_file,
        "w",
        newline=""
    ) as f:

        writer = csv.writer(f)

        writer.writerow([
            "encoding",
            "configuration",
            "dna_cost_nt_per_bit",
            "error_tolerance",
            "success_rate_at_tolerance"
        ])

        for encoding, config_data in cost_vs_tolerance.items():

            for config_name, result in config_data.items():

                writer.writerow([
                    encoding,
                    config_name,
                    result["dna_cost"],
                    result["error_tolerance"],
                    result["success_rate"]
                ])

    print(f"Saved summary: {summary_file}")
    # ================================================================
    # PRINT COST / ROBUSTNESS RESULTS
    # ================================================================

    print("\n")
    print("=" * 80)
    print(
        f"DNA COST VS ERROR TOLERANCE "
        f"(target recovery = {target_recovery:.0%})"
    )
    print("=" * 80)


    for encoding, config_data in cost_vs_tolerance.items():

        print(f"\n{encoding}")

        for config_name, result in config_data.items():

            dna_cost = result['dna_cost']
            tolerance = result['error_tolerance']
            success_rate = result['success_rate']

            if tolerance is None:

                tolerance_text = (
                    "not reached within tested range"
                )

                success_text = "-"

            else:

                tolerance_text = (
                    f"{tolerance:.5f}"
                )

                success_text = (
                    f"{success_rate:.1%}"
                )

            print(
                f"  {config_name:<45} "
                f"cost={dna_cost:.4f} nt/bit   "
                f"tolerance={tolerance_text:<35} "
                f"success={success_text}"
            )


    # ================================================================
    # CREATE SCATTER PLOT
    # ================================================================

    plt.figure(figsize=(12, 8))


    # Short names for the different encodings
    encoding_prefix = {
        'goldman': 'G',
        'church': 'C',
        'wukong': 'W',
        'gcplus': 'GC',
        'hedges': 'H',
        'max_density': 'MD',
        'no_homopolymer': 'NH'
    }


    # ================================================================
    # LOOP THROUGH ALL ENCODINGS
    # ================================================================

    for encoding, config_data in cost_vs_tolerance.items():

        x = []
        y = []
        labels = []

        prefix = encoding_prefix.get(
            encoding,
            encoding
        )

        # ------------------------------------------------------------
        # LOOP THROUGH ALL CONFIGURATIONS OF THIS ENCODING
        # ------------------------------------------------------------

        for config_number, (config_name, result) in enumerate(
            config_data.items()
        ):

            tolerance = result['error_tolerance']

            # --------------------------------------------------------
            # Configurations that never reach the target recovery
            # cannot be plotted at a meaningful tolerance value.
            # --------------------------------------------------------

            if tolerance is None:
                continue

            x.append(
                result['dna_cost']
            )

            y.append(
                tolerance
            )

            # --------------------------------------------------------
            # Short label for this configuration
            #
            # Examples:
            # G1, G2, G3
            # C1, C2, C3
            # W1, W2, W3
            # GC1, GC2, GC3
            # --------------------------------------------------------

            if encoding in ['church', 'wukong']:

                redundancy = int(
                    result['configuration']['add_redundancy']
                )

                rs_num = result['configuration']['rs_num']

                labels.append(
                    f"{prefix}{redundancy}-{rs_num}"
                )

            else:

                labels.append(
                    f"{prefix}{config_number + 1}"
                )


        # ------------------------------------------------------------
        # PLOT ALL CONFIGURATIONS FOR THIS ENCODING
        # ------------------------------------------------------------

        if len(x) > 0:

            plt.scatter(
                x,
                y,
                s=80,
                label=encoding
            )

            # --------------------------------------------------------
            # Add short labels next to each dot
            # --------------------------------------------------------
            label_offsets = [
                (7, 7),
                (7, -12),
                (-7, 7),
                (-7, -12),
                (12, 0),
                (-12, 0),
                (0, 12),
                (0, -15),
            ]

            for label_number, (x_value, y_value, label) in enumerate(
                zip(x, y, labels)
            ):

                offset = label_offsets[label_number % len(label_offsets)]

                plt.annotate(
                    label,
                    (x_value, y_value),
                    xytext=offset,
                    textcoords='offset points',
                    fontsize=8,
                    ha='center',
                    va='center',
                    bbox=dict(
                        boxstyle='round,pad=0.2',
                        facecolor='white',
                        edgecolor='none',
                        alpha=0.8
                    )
                )


    # ================================================================
    # AXES AND TITLE
    # ================================================================

    plt.xlabel(
        "DNA storage cost (nt/bit)"
    )

    plt.ylabel(
        f"Maximum tested error rate with ≥ "
        f"{target_recovery:.0%} recovery"
    )

    plt.title(
        "DNA storage cost vs sequencing-error tolerance"
    )

    plt.grid(
        True,
        alpha=0.3
    )

    plt.legend()

    plt.tight_layout()


    # ================================================================
    # SAVE PNG
    # ================================================================

    cost_png = os.path.join(
        output_dir,
        "dna_cost_vs_error_tolerance.png"
    )

    plt.savefig(
        cost_png,
        dpi=300,
        bbox_inches='tight'
    )

    print(
        f"Saved: {cost_png}"
    )


    # ================================================================
    # SAVE PDF
    # ================================================================

    cost_pdf = os.path.join(
        output_dir,
        "dna_cost_vs_error_tolerance.pdf"
    )

    plt.savefig(
        cost_pdf,
        bbox_inches='tight'
    )

    print(
        f"Saved: {cost_pdf}"
    )


    plt.show()

    plt.close()
# if __name__ == '__main__':
#     PARAMETER_SWEEPS = {


#     'church': {
#         'add_redundancy': [False, True],
#         'rs_num': [0, 1, 2, 3, 4, 5],
#     },

#     'wukong': {
#         'add_redundancy': [False, True],
#         'rs_num': [0, 1, 2, 3, 4, 5],
#     },

#     'gcplus': {
#         'gcplus_c1': [1, 2, 3, 4, 5, 6],
#     },

#     'hedges': {
#         'hedges_coderate': [1, 2, 3, 4, 5, 6],
#     },

#     'max_density': {
#         'reed_solo_percentage': [
#             0.50, 0.60, 0.70, 0.80, 0.90, 1.00
#         ],
#     },

#     'no_homopolymer': {
#         'reed_solo_percentage': [
#             0.50, 0.60, 0.70, 0.80, 0.90, 1.00
#         ],
#     },
# }
    
#     encodings = ['goldman', 'church', 'no_homopolymer']

#     base_params = {
#         'filename': 'testfilesimsall.txt',
#         'binarization_method': 'default',
#         'synthesis_method': 'nosynthpoly',
#         'sequencing_method': 'iid',
#         'recovery_method': 'pass_through',
#         'min_coverage': 1,
#         'clustering_method': 'pass_through',
#         'storage_conditions': None,
#         'kmer_seed': 42,
#         'similarity_threshold': 1,
#         'gap_penalty': 1.0,
#     }

#     no_error_params = {
#         'mean': 1,
#         'std_dev': 0,
#         "iid_error_rate": 0.0,
#         "iid_substitution_rate": 0.0,
#         "iid_insertion_rate": 0.0,
#         "iid_deletion_rate": 0.0,
#     }
    
#     error_rates = [0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01]  # Example error rates to test

#     defoultparams = {}      
#     for encoding in encodings:
#         parameters=base_params.copy()
#         parameters.update(default_values_of_encoding(encoding))
#         parameters.update(no_error_params)
#         if encoding == 'no_homopolymer' or encoding == 'max_density':
#             parameters['clustering_method'] = None
#             parameters['recovery_method'] = None
#         print(f"Parameters for {encoding}: {parameters}")
#         defoultparams[encoding] = Params(**parameters)
        

#     # noECencodinglengths = encode_and_count(encodings, defoultparams)

#     # print(noECencodinglengths)

#     # final_params = encode_and_count_with_error_correction(encodings, defoultparams, noECencodinglengths)

#     noECencodinglengths = encode_and_count(encodings, defoultparams)

#     print(noECencodinglengths)

#     final_params = encode_and_count_with_error_correction(
#         encodings,
#         defoultparams,
#         noECencodinglengths
# )


#     encodings = ['goldman', 'church', 'no_homopolymer']

#     sim_resoults = simulate_cost_analysis(encodings, final_params, error_rates)

#     print("Simulation results for different error rates:")
#     for encoding, error_data in sim_resoults.items():

#         print(f"\nEncoding: {encoding}")

#         for error_rate, runs in error_data.items():

#             success_count = sum(
#                 run['success']
#                 for run in runs
#             )

#             success_rate = (
#                 success_count / len(runs)
#             )

#             dna_length = runs[0]['dna_length']
#             information_bits = runs[0]['information_bits']
#             dna_cost = runs[0]['dna_cost_nt_per_bit']
#             expected_errors = runs[0]['expected_errors']

#             print(
#                 f"  Error Rate: {error_rate:.4g}, "
#                 f"Success: {success_count}/{len(runs)} "
#                 f"({success_rate:.1%}), "
#                 f"DNA: {dna_length} nt, "
#                 f"Cost: {dna_cost:.3f} nt/bit, "
#                 f"Expected errors: {expected_errors:.2f}"
#             )

#     # Save the results to a file for further analysis
#     output_file = "simulation_results.txt"
#     with open(output_file, "w") as f:
#         for encoding, error_data in sim_resoults.items():
#             f.write(f"Encoding: {encoding}\n")
#             for error_rate, runs in error_data.items():
#                 success_count = sum(run['success'] for run in runs)
#                 f.write(f"  Error Rate: {error_rate}, Successful Runs: {success_count}/{len(runs)}\n")

#     print(f"Simulation results saved to {output_file}")

#     #make plot with the results
#     import matplotlib.pyplot as plt

#     for encoding, error_data in sim_resoults.items():
#         error_rates = sorted(error_data.keys())
#         success_rates = [sum(run['success'] for run in error_data[er]) / len(error_data[er]) for er in error_rates]
#         plt.plot(error_rates, success_rates, 'o-', label=encoding)

#     plt.xlabel('Error Rate')
#     plt.ylabel('Success Rate')
#     plt.title('Simulation Results for Different Error Rates with nearly equal bp count for all encodings')
#     plt.legend()
#     plt.savefig('simulation_results_plot.png')
#     plt.show()


#     # ---------------------------------------------------------
#     # COST VS ERROR TOLERANCE
#     # ---------------------------------------------------------

#     target_recovery = 0.90

#     cost_vs_tolerance = {}

#     for encoding, error_data in sim_resoults.items():

#         # Get DNA cost from the first available simulation
#         first_error_rate = next(iter(error_data))
#         first_run = error_data[first_error_rate][0]

#         dna_cost = first_run["dna_cost_nt_per_bit"]

#         # Sort error rates from low -> high
#         sorted_error_rates = sorted(error_data.keys())

#         tolerated_error_rate = None
#         tolerated_success_rate = None

#         for error_rate in sorted_error_rates:

#             runs = error_data[error_rate]

#             success_count = sum(
#                 run["success"]
#                 for run in runs
#             )

#             success_rate = success_count / len(runs)

#             # Keep going while we satisfy the target
#             if success_rate >= target_recovery:
#                 tolerated_error_rate = error_rate
#                 tolerated_success_rate = success_rate
#             else:
#                 # Since error rates are increasing,
#                 # we can stop here.
#                 break

#         cost_vs_tolerance[encoding] = {
#             "dna_cost": dna_cost,
#             "error_tolerance": tolerated_error_rate,
#             "success_rate": tolerated_success_rate
#         }


#     # Print results
#     print("\nCost vs. error tolerance:")

#     for encoding, result in cost_vs_tolerance.items():

#         print(
#             f"{encoding}: "
#             f"{result['dna_cost']:.3f} nt/bit, "
#             f"max error rate = "
#             f"{result['error_tolerance']}, "
#             f"success rate = {result['success_rate']:.1%}"
#         )


#     x = []
#     y = []
#     labels = []

#     for encoding, result in cost_vs_tolerance.items():

#         if result["error_tolerance"] is None:
#             continue

#         x.append(result["dna_cost"])
#         y.append(result["error_tolerance"])
#         labels.append(encoding)

#     plt.figure(figsize=(8, 6))

#     plt.scatter(x, y)

#     for i, label in enumerate(labels):
#         plt.annotate(
#             label,
#             (x[i], y[i]),
#             xytext=(5, 5),
#             textcoords="offset points"
#         )

#     plt.xlabel("DNA cost (nt/bit)")
#     plt.ylabel("Maximum tolerated error rate")
#     plt.title(
#         f"DNA cost vs. error tolerance "
#         f"({target_recovery:.0%} recovery threshold)"
#     )

#     plt.grid(True)
#     plt.tight_layout()
#     plt.savefig("cost_vs_error_tolerance.png", dpi=300)
#     plt.show()