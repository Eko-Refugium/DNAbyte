
import os

from dnabyte.params import Params
from dnabyte.binarize import Binarize
from simulations.simulation import Simulation
from dnabyte.data_classes.base import Data
from dnabyte.encode import Encode

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
    target_bp = max(noECencodinglengths.values())
    max_encoding = max(noECencodinglengths, key=noECencodinglengths.get)
    final_encodedbp[max_encoding] = target_bp
    final_encodedbpnomean[max_encoding] = target_bp
    encodings.remove(max_encoding)
    for encoding in encodings:
        if encoding == 'goldman' or encoding == 'yinyang':
            if encoding == 'yinyang':
                defultparams[encoding].mean = 5
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding failed for {encoding}")
            while True:
                if success < target_bp+300:
                    defultparams[encoding].mean += 1
                elif success > target_bp+300:
                    defultparams[encoding].mean -= 1
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean
            
        elif encoding == 'church' or encoding == 'wukong':
            defultparams[encoding].add_redundancy = True
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding add_redundancy failed for {encoding}")
            if success > target_bp+300:
                defultparams[encoding].add_redundancy = False
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding add_redundancy failed for {encoding}")
            while defultparams[encoding].rs_num < 10:
                if success < target_bp:
                    defultparams[encoding].rs_num += 1
                elif success > target_bp:
                    defultparams[encoding].rs_num -= 1
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding add_redundancy failed for {encoding}")

            while True:
                if success < target_bp+300:
                    defultparams[encoding].mean += 1
                elif success > target_bp+300:
                    defultparams[encoding].mean -= 1
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean

        elif encoding == 'gcplus':
            defultparams[encoding].gcplus_c1 = 1
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding gcplus_c1 failed for {encoding}")
            while defultparams[encoding].gcplus_c1 < 10:
                if success < target_bp+300:
                    defultparams[encoding].gcplus_c1 += 1
                elif success > target_bp+300:
                    defultparams[encoding].gcplus_c1 -= 1
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding gcplus_c1 failed for {encoding}")
            while True:
                if success < target_bp+300:
                    defultparams[encoding].mean += 1
                elif success > target_bp+300:
                    defultparams[encoding].mean -= 1
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding failed for {encoding}")
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
            while (defultparams[encoding].reed_solo_percentage > 0.5):
                if success < target_bp+300:
                    defultparams[encoding].reed_solo_percentage = max(defultparams[encoding].reed_solo_percentage - 0.01, 0.0)
                elif success > target_bp+300:
                    defultparams[encoding].reed_solo_percentage = min(defultparams[encoding].reed_solo_percentage + 0.01, 1.0)
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding failed for {encoding}")
            while True:
                if success < target_bp+300:
                    defultparams[encoding].mean += 1
                elif success > target_bp+300:
                    defultparams[encoding].mean -= 1
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding failed for {encoding}")
            final_encodedbp[encoding] = success
            final_encodedbpnomean[encoding] = encodedbpnomean

        elif encoding == 'hedeges':
            defultparams[encoding].hedges_coderate = 3
            success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
            if success is None:
                raise RuntimeError(f"Encoding hedges failed for {encoding}")
            while defultparams[encoding].hedges_coderate < 10:
                if success < target_bp+300:
                    defultparams[encoding].hedges_coderate += 1
                elif success > target_bp+300:
                    defultparams[encoding].hedges_coderate -= 1
                    break
                success, encodedbpnomean = encode_and_count_single_encoding(encoding, defultparams)
                if success is None:
                    raise RuntimeError(f"Encoding failed for {encoding}")
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
            defultparams[encoding].iid_insertion_rate = error_rate
            defultparams[encoding].iid_deletion_rate = error_rate
            for i in range(3):  # Run each simulation 20 times for averaging
                print(f"Run {i+1} for encoding {encoding} at error rate {error_rate}")
                sim = Simulation([defultparams[encoding]])
                results = sim.run()
                
                # Check if any simulation succeeded
                encsuccess = any(
                    sim_result.get('status') == 'SUCCESS' 
                    for sim_result in results.values()
                )
                
                if encoding not in error_simulation_results:
                    error_simulation_results[encoding] = {}
                
                if error_rate not in error_simulation_results[encoding]:
                    error_simulation_results[encoding][error_rate] = []
                
                error_simulation_results[encoding][error_rate].append({
                    "run": i + 1,
                    "success": encsuccess,
                    "results": results
                })
            
            # Check if any simulation succeeded
            encsuccess = any(
                sim_result.get('status') == 'SUCCESS' 
                for sim_result in results.values()
            )

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


if __name__ == '__main__':
    encodings = ['goldman', 'church', 'gcplus', 'max_density', 'no_homopolymer', 'wukong', 'yinyang', 'hedges']

    base_params = {
        'filename': 'Bohemian_Rhapsody_Lyrics.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'iid',
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
        "iid_error_rate": 0.0,
        "iid_substitution_rate": 0.0,
        "iid_insertion_rate": 0.0,
        "iid_deletion_rate": 0.0,
    }
    
    error_rates = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.12, 0.16, 0.20]

    defoultparams = {}      
    for encoding in encodings:
        parameters=base_params.copy()
        parameters.update(default_values_of_encoding(encoding))
        parameters.update(no_error_params)
        print(f"Parameters for {encoding}: {parameters}")
        defoultparams[encoding] = Params(**parameters)
        

    noECencodinglengths = encode_and_count(encodings, defoultparams)

    final_params = encode_and_count_with_error_correction(encodings, defoultparams, noECencodinglengths)

    sim_resoults = simulate_cost_analysis(encodings, final_params, error_rates)

    print("Simulation results for different error rates:")
    for encoding, error_data in sim_resoults.items():
        print(f"\nEncoding: {encoding}")
        for error_rate, runs in error_data.items():
            success_count = sum(run['success'] for run in runs)
            print(f"  Error Rate: {error_rate}, Successful Runs: {success_count}/{len(runs)}")


    
    