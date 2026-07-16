"""
Generic encoding analyzer for calculating encoding-specific properties.
Produces metrics that can be used to tune simulation parameters.
"""
import math

import numpy as np
from dnabyte.params import Params
from dnabyte.data_classes.base import Data
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode


class AnalyzeEncoding:
    """
    Analyze encoding properties to support parameter tuning in simulations.
    
    Calculates metrics for a given encoding that can be used to:
    - Set target copy numbers (mean/std_dev)
    - Predict synthesis requirements
    - Compare encoding efficiency
    - Plan error correction strategies
    """
    
    def __init__(self, params):
        """
        Initialize analyzer with encoding parameters.
        
        :param params: Params object with encoding_method, codeword_length, etc.
        """
        self.params = params
        self.metrics = {}
    
    def calculate_data_density(self):
        """
        Calculate data density: bits of information per DNA base.
        
        :return: bits per base
        """
        # Create dummy input data (1000 bits)
        dummy_data = '0' * 1000
        
        # Binarize
        bin_obj = Binarize(self.params)
        data_obj = Data(file_paths=[])
        binary_code = type('obj', (object,), {'data': dummy_data, 'data_type': 'binary'})()
        
        # Encode
        coder = Encode(self.params)
        data_enc, enc_info = coder.encode(binary_code)
        
        # Calculate: bits input / bases output
        total_bases = sum(len(seq) for seq in data_enc.data)
        bits_input = len(dummy_data)
        
        density = bits_input / total_bases if total_bases > 0 else 0
        self.metrics['data_density'] = density
        
        return density
    
    def calculate_synthesis_cost(self):
        from dnabyte.encoding.gcplus.src.GCPdna.GCP_Encode_DNA import GCP_Encode_DNA_brute
        from dnabyte.encoding.gcplus.encode import _load_codebook
        params = self.params
        if params.encoding_method == "church":
            dna_len = params.codeword_length
            if params.add_primer:
                dna_len -= 2 * params.primer_length

            # Step 2: binary capacity
            capacity = math.ceil(dna_len * params.bit_per_base)

            
            bytes_num = math.ceil(capacity / 8)
            rs_group = math.ceil(bytes_num / 255)
            # Step 3: reserve RS bits
            capacity -= params.rs_num * 8 * rs_group

            # Initial guess
            index_len = 1

            bin = Binarize(params)

            # Handle both file_paths (for compressed) and filename (for default binarization)
            # Prepend ./tests/testfiles/ like testbase does
            if hasattr(params, 'file_paths') and params.file_paths:
                file_paths = ['./tests/testfiles/' + fp for fp in params.file_paths]
            elif hasattr(params, 'filename'):
                file_paths = ['./tests/testfiles/' + params.filename]
            else:
                raise ValueError("params must have either 'file_paths' or 'filename' attribute")
            
            data_obj = Data(file_paths=file_paths)
            binary_code = bin.binarize(data_obj)

            while True:

                payload = capacity - index_len

                if payload <= 0:
                    raise ValueError("Payload <= 0")

                # Data strands only
                

                data_segments = math.ceil(binary_code.data_bits / payload)

                # Apply Storage-D redundancy
                if params.add_redundancy:
                    total_segments = data_segments + data_segments // 2
                else:
                    total_segments = data_segments

                new_index = math.ceil(math.log2(total_segments))

                if new_index == index_len:
                    break

                index_len = new_index

            return {
                "payload_bits": payload,
                "index_bits": index_len,
                "data_segments": data_segments,
                "total_strands": total_segments
            }
        
        if params.encoding_method == "goldman":
            bin = Binarize(params)

            # Handle both file_paths (for compressed) and filename (for default binarization)
            # Prepend ./tests/testfiles/ like testbase does
            if hasattr(params, 'file_paths') and params.file_paths:
                file_paths = ['./tests/testfiles/' + fp for fp in params.file_paths]
            elif hasattr(params, 'filename'):
                file_paths = ['./tests/testfiles/' + params.filename]
            else:
                raise ValueError("params must have either 'file_paths' or 'filename' attribute")
            
            data_obj = Data(file_paths=file_paths)
            binary_code = bin.binarize(data_obj)

            total_bits = binary_code.data_bits
            ternary_length = math.ceil(total_bits / 8) * 5

            index_len = 1
            

            while True:

                # payload after reserving ternary index
                payload = params.codeword_length - index_len

                if payload <= 0:
                    raise ValueError("Sequence too short.")

                # Goldman slides by 1/4 of the payload
                step = payload // 4

                if step == 0:
                    raise ValueError("Payload too small.")

                # number of overlapping segments
                if ternary_length <= payload:
                    segments = 1
                else:
                    segments = math.ceil((ternary_length - payload) / step) + 1

                new_index = math.ceil(math.log(segments, 3))

                if new_index == index_len:
                    break

                index_len = new_index

            return {
                "ternary_symbols": ternary_length,
                "payload_per_strand": payload,
                "step": step,
                "index_length": index_len,
                "segments": segments
            }
        if params.encoding_method == "wukong":
            dna_len = params.codeword_length

            if params.add_primer:
                dna_len -= 2 * params.primer_length

            # Wukong stores 2 bits/base
            capacity = math.ceil(dna_len * params.bit_per_base)

            # Reserve RS parity
            bytes_num = math.ceil(capacity / 8)
            rs_group = math.ceil(bytes_num / 255)
            capacity -= params.rs_num * 8 * rs_group

            # Wukong reserves one extra index bit
            index_redundancy = 1
            index_len = 1 + index_redundancy

            while True:

                payload = capacity - index_len

                if payload <= 0:
                    raise ValueError("Payload <= 0")

                # Number of binary segments
                data_segments = math.ceil(binary_code.data_bits / payload)

                # Optional Storage-D redundancy
                if params.add_redundancy:
                    total_segments = data_segments + data_segments // 2
                else:
                    total_segments = data_segments

                # Index must identify every binary segment
                new_index = math.ceil(math.log2(total_segments)) + index_redundancy

                if new_index == index_len:
                    break

                index_len = new_index

            # Wukong packs two binary segments into one DNA strand
            dna_strands = math.ceil(total_segments / 2)

            return {
                "payload_bits": payload,
                "index_bits": index_len,
                "data_segments": data_segments,
                "total_segments": total_segments,
                "dna_strands": dna_strands,
            }
        
        if params.encoding_method == "gcplus":

            bin = Binarize(params)

            if hasattr(params, 'file_paths') and params.file_paths:
                file_paths = ['./tests/testfiles/' + fp for fp in params.file_paths]
            elif hasattr(params, 'filename'):
                file_paths = ['./tests/testfiles/' + params.filename]
            else:
                raise ValueError(
                    "params must have either file_paths or filename attribute"
                )

            data_obj = Data(file_paths=file_paths)
            binary_code = bin.binarize(data_obj)

            total_bits = binary_code.data_bits


            # GC+ parameters
            k = int(getattr(params, "gcplus_k", 168))
            l = int(getattr(params, "gcplus_l", 8))
            c1 = int(getattr(params, "gcplus_c1", 2))


            # Calculate oligo count
            codewords = math.ceil(total_bits / k)


            # Calculate actual GC+ DNA length
            codebook = _load_codebook()

            dummy = [0] * k

            dna, dna_len, N, K, q, U, X, check_par = GCP_Encode_DNA_brute(
                dummy,
                l,
                c1,
                codebook
            )


            if dna_len > params.sequence_length:
                raise ValueError(
                    f"GC+ oligo length {dna_len} exceeds "
                    f"requested sequence_length {params.sequence_length}"
                )


            return {
                "payload_bits": k,
                "sequence_length": dna_len,
                "codewords": codewords,
                "total_bits": total_bits,

                # GC+ internal parameters
                "K": K,
                "N": N,
                "q": q,

                "gcplus_k": k,
                "gcplus_l": l,
                "gcplus_c1": c1
            }

            
    def calculate_coverage_requirement(self, target_reliability=0.99):
        """
        Calculate recommended oligonucleotide copy number for target reliability.
        
        Based on Poisson distribution of copy coverage.
        
        :param target_reliability: Probability (0-1) of at least one copy being error-free
        :return: Recommended mean copy number
        """
        from scipy.special import gammainc, gamma
        
        # Poisson parameter λ for getting P(X >= 1) = target_reliability
        # P(X >= 1) = 1 - P(X = 0) = 1 - e^(-λ)
        # Solving: e^(-λ) = 1 - target_reliability
        
        lambda_mean = -np.log(1 - target_reliability)
        
        self.metrics['coverage_mean'] = lambda_mean
        self.metrics['target_reliability'] = target_reliability
        
        return lambda_mean
    
    def calculate_error_correction_capability(self):
        """
        Calculate error correction capability based on encoding AND simulation parameters.
        
        NOTE: Error correction capability depends on:
        - Encoding method (base resilience)
        - Mean/std_dev (copy number redundancy)
        - Synthesis method (MESA increases copies, assembly is complex)
        - Storage duration (degradation reduces effective coverage)
        - Assembly probability (affects successful recovery)
        
        This method provides BASE ESTIMATES - should be MEASURED via benchmarks.
        
        :return: Estimated maximum tolerable error rate (0-1)
        """
        encoding_method = self.params.encoding_method
        
        # Base encoding resilience (WITHOUT considering other parameters)
        base_resilience = {
            'max_density': 0.15,      # Most robust DNA encoding
            'church': 0.10,           # Good error correction
            'wukong': 0.05,           # Limited error correction
            'linear_chain': 0.12,
            'poly_chain': 0.08,
            'linear_binom': 0.10,
            'poly_binom': 0.07,
        }
        
        max_error_base = base_resilience.get(encoding_method, 0.10)
        
        # Adjust for synthesis method (increases/decreases redundancy)
        synthesis_adjustment = 1.0
        synthesis_method = getattr(self.params, 'synthesis_method', None)
        if synthesis_method == 'mesa':
            synthesis_adjustment = 1.3  # MESA creates multiple copies
        elif synthesis_method == 'assembly':
            synthesis_adjustment = 0.9  # Assembly reduces effective coverage
        
        # Adjust for copy numbers (redundancy protection)
        mean_copies = getattr(self.params, 'mean', 10)
        if mean_copies >= 20:
            copy_adjustment = 1.2
        elif mean_copies >= 10:
            copy_adjustment = 1.0
        else:
            copy_adjustment = 0.8
        
        # Adjust for storage duration (degradation)
        years = getattr(self.params, 'years', 10)
        if years <= 1:
            storage_adjustment = 1.0
        elif years <= 100:
            storage_adjustment = 0.95
        elif years <= 1000:
            storage_adjustment = 0.8
        else:
            storage_adjustment = 0.6
        
        # Adjust for assembly probability (affects recovery success)
        assembly_prob = getattr(self.params, 'assembly_probability', 1.0)
        if assembly_prob < 0.9:
            assembly_adjustment = 0.85
        elif assembly_prob < 0.95:
            assembly_adjustment = 0.95
        else:
            assembly_adjustment = 1.0
        
        # Combined estimate
        max_error_rate = max_error_base * synthesis_adjustment * copy_adjustment * storage_adjustment * assembly_adjustment
        
        # Store breakdown for transparency
        self.metrics['max_error_rate'] = max_error_rate
        self.metrics['error_correction_breakdown'] = {
            'base_resilience': max_error_base,
            'synthesis_adjustment': synthesis_adjustment,
            'copy_adjustment': copy_adjustment,
            'storage_adjustment': storage_adjustment,
            'assembly_adjustment': assembly_adjustment,
        }
        
        return max_error_rate
    
    def calculate_assembly_efficiency(self):
        """
        Estimate assembly efficiency needed based on encoding AND synthesis parameters.
        
        NOTE: Assembly probability requirement depends on:
        - Encoding complexity (simpler → lower requirement)
        - Synthesis method (assembly method more sensitive to probability)
        - Number of synthesis steps needed
        - Whether there are alternative recovery paths
        
        This provides ESTIMATES - actual requirements should be MEASURED via benchmarks.
        
        :return: Recommended assembly probability (0-1)
        """
        encoding_method = self.params.encoding_method
        synthesis_method = getattr(self.params, 'synthesis_method', None)
        
        # Base assembly requirement by encoding complexity
        base_requirements = {
            'max_density': 0.85,      # Simple encoding
            'church': 0.90,           # Moderate
            'wukong': 0.95,           # Complex
            'linear_chain': 0.80,     # Simpler assembly
            'poly_chain': 0.75,       # More permissive
            'linear_binom': 0.82,
            'poly_binom': 0.78,
        }
        
        required_assembly_prob = base_requirements.get(encoding_method, 0.90)
        
        # Adjust based on synthesis method
        if synthesis_method == 'assembly':
            required_assembly_prob += 0.05  # Assembly critical for success
        elif synthesis_method == 'mesa':
            required_assembly_prob -= 0.05  # MESA provides redundancy
        
        # Clamp to [0, 1]
        required_assembly_prob = max(0.0, min(1.0, required_assembly_prob))
        
        self.metrics['required_assembly_probability'] = required_assembly_prob
        self.metrics['assembly_efficiency_notes'] = {
            'encoding': encoding_method,
            'synthesis_method': synthesis_method or 'none',
            'base_requirement': base_requirements.get(encoding_method, 0.90),
        }
        
        return required_assembly_prob
    
    def calculate_all_metrics(self):
        """
        Calculate all available metrics for this encoding.
        
        :return: Dictionary of metrics
        """
        try:
            self.calculate_data_density()
        except:
            self.metrics['data_density'] = None
        
        self.calculate_synthesis_cost()
        self.calculate_coverage_requirement()
        self.calculate_error_correction_capability()
        self.calculate_assembly_efficiency()
        
        return self.metrics
    
    def print_report(self):
        """Print a human-readable report of all metrics with parameter dependencies."""
        if not self.metrics:
            self.calculate_all_metrics()
        
        print(f"\n{'='*70}")
        print(f"ENCODING ANALYSIS REPORT: {self.params.encoding_method}")
        print(f"{'='*70}\n")
        
        # Parameter context
        print("SIMULATION PARAMETERS:")
        print(f"  Synthesis method: {getattr(self.params, 'synthesis_method', 'none')}")
        print(f"  Mean copies: {getattr(self.params, 'mean', 10)}")
        print(f"  Storage duration: {getattr(self.params, 'years', 10)} years")
        print(f"  Assembly probability: {getattr(self.params, 'assembly_probability', 1.0):.2f}")
        print()
        
        print(f"Data Density:")
        if self.metrics.get('data_density'):
            print(f"  {self.metrics['data_density']:.4f} bits/base")
        else:
            print(f"  N/A (requires full encoding test)")
        
        print(f"\nSynthesis Cost (for synthesis-based approaches):")
        if self.metrics.get('synthesis_cost') == 0:
            print(f"  {self.metrics.get('synthesis_cost_notes', 'N/A (no synthesis)')}")
        elif self.metrics.get('synthesis_cost') is not None:
            print(f"  {self.metrics.get('synthesis_cost', 'N/A'):.0f} bases to synthesize per Kbit input data")
            if self.metrics.get('synthesis_cost_breakdown'):
                bd = self.metrics['synthesis_cost_breakdown']
                print(f"  Breakdown (test encoding of {bd['input_bits']} bits):")
                print(f"    Codewords created: {bd['num_codewords']}")
                print(f"    Total bases: {bd['total_bases']}")
                print(f"    Bases per Kbit input: {bd['bases_per_kbit']:.0f}")
            else:
                print(f"  {self.metrics.get('synthesis_cost_notes', 'Estimated from parameters')}")
        else:
            print(f"  {self.metrics.get('synthesis_cost_notes', 'N/A')}")

        
        print(f"\nCoverage Requirement (for {self.metrics.get('target_reliability', 0.99):.0%} reliability):")
        print(f"  Mean copy number: {self.metrics.get('coverage_mean', 'N/A'):.2f}")
        
        print(f"\nError Correction Capability:")
        print(f"  Max tolerable error rate: {self.metrics.get('max_error_rate', 'N/A'):.1%}")
        if self.metrics.get('error_correction_breakdown'):
            bd = self.metrics['error_correction_breakdown']
            print(f"  Breakdown:")
            print(f"    Base resilience: {bd['base_resilience']:.2f}")
            print(f"    × Synthesis adjustment: {bd['synthesis_adjustment']:.2f}")
            print(f"    × Copy redundancy: {bd['copy_adjustment']:.2f}")
            print(f"    × Storage degradation: {bd['storage_adjustment']:.2f}")
            print(f"    × Assembly reliability: {bd['assembly_adjustment']:.2f}")
            print(f"    = Final estimate: {self.metrics.get('max_error_rate', 'N/A'):.1%}")
        
        print(f"\nAssembly Requirements:")
        print(f"  Required assembly probability: {self.metrics.get('required_assembly_probability', 'N/A'):.2f}")
        if self.metrics.get('assembly_efficiency_notes'):
            notes = self.metrics['assembly_efficiency_notes']
            print(f"  Context: {notes['synthesis_method']} synthesis for {notes['encoding']}")
        
        print(f"\n{'='*70}")
        print("NOTE: These are ESTIMATES based on parameters. Validate with benchmarks!")
        print(f"{'='*70}\n")

    
    def recommend_parameters(self):
        """
        Recommend simulation parameters based on calculated metrics.
        
        :return: Dictionary of recommended parameters
        """
        if not self.metrics:
            self.calculate_all_metrics()
        
        recommendations = {
            'assembly_probability': self.metrics.get('required_assembly_probability', 0.90),
            'mean_copies': int(np.ceil(self.metrics.get('coverage_mean', 10))),
            'max_sequencing_error_rate': self.metrics.get('max_error_rate', 0.10),
        }
        
        return recommendations
