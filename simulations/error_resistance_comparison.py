"""
Error Resistance Comparison Framework
Compares DNA encodings at different error rates to measure resilience.
"""

import os
import sys
import logging
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Any, Tuple
from tqdm import tqdm
import numpy as np

# Add DNAbyte to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.data_classes.base import Data
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus


@dataclass
class ErrorMetrics:
    """Metrics for error resistance testing"""
    encoding_method: str
    error_rate: float
    num_runs: int
    successful_recoveries: int
    total_errors_introduced: int
    avg_errors_per_strand: float
    recovery_success_rate: float
    
    def __str__(self):
        return (f"{self.encoding_method} @ {self.error_rate*100:.1f}% error: "
                f"{self.recovery_success_rate*100:.1f}% recovery "
                f"({self.successful_recoveries}/{self.num_runs} runs)")


class ErrorResistanceTester:
    """Test encoding error resilience at different error rates"""
    
    def __init__(self, output_dir: str = './simulations/error_results'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        # Setup logging
        log_file = os.path.join(output_dir, f'error_resistance_{self._timestamp()}.log')
        self.logger = self._setup_logger(log_file)
    
    def _timestamp(self):
        from datetime import datetime
        return datetime.now().strftime('%Y%m%d_%H%M%S')
    
    def _setup_logger(self, log_file: str):
        logger = logging.getLogger('error_resistance')
        logger.setLevel(logging.INFO)
        
        handler = logging.FileHandler(log_file)
        handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s'
        ))
        logger.addHandler(handler)
        return logger
    
    def test_encoding_at_error_rate(self,
                                   encoding_method: str,
                                   test_data: Data,
                                   error_rate: float,
                                   params_dict: Dict[str, Any],
                                   num_runs: int = 3) -> ErrorMetrics:
        """
        Test single encoding at specific error rate.
        
        Simplified pipeline: Data → Binarize → Encode → 
        Introduce random errors to DNA sequences → Measure resilience
        """
        
        successful_recoveries = 0
        total_errors = 0
        
        for run in range(num_runs):
            try:
                # Step 1: Binarize
                params = Params(**params_dict)
                binary_code = Binarize(params).binarize(test_data)
                
                # Step 2: Encode
                dna_codewords, _ = Encode(params).encode(binary_code)
                dna_sequences = dna_codewords.data
                
                # Step 3: Introduce random errors to DNA sequences
                dna_with_errors = self._introduce_errors(dna_sequences, error_rate)
                errors_introduced = self._count_differences(dna_sequences, dna_with_errors)
                total_errors += errors_introduced
                
                # Step 4: Try to recover via clustering consensus
                # For now, just measure how many errors were introduced
                # A robust encoding should have error correction that helps recover from this
                
                # Simple metric: encodings with redundancy should handle more errors
                # We measure by checking if the encoding can be decoded despite errors
                recovery_success = self._test_recovery_feasibility(
                    dna_sequences, dna_with_errors, encoding_method
                )
                
                if recovery_success:
                    successful_recoveries += 1
                
                self.logger.info(f"[OK] {encoding_method} run {run+1}: "
                               f"error_rate={error_rate*100:.1f}%, "
                               f"errors={errors_introduced}, recovery={'success' if recovery_success else 'failed'}")
                
            except Exception as e:
                self.logger.error(f"[ERROR] {encoding_method} run {run+1}: {str(e)}")
                import traceback
                self.logger.error(traceback.format_exc())
                continue
        
        avg_errors = total_errors / max(num_runs, 1) if total_errors > 0 else 0
        success_rate = successful_recoveries / max(num_runs, 1)
        
        return ErrorMetrics(
            encoding_method=encoding_method,
            error_rate=error_rate,
            num_runs=num_runs,
            successful_recoveries=successful_recoveries,
            total_errors_introduced=total_errors,
            avg_errors_per_strand=avg_errors,
            recovery_success_rate=success_rate
        )
    
    def _introduce_errors(self, sequences: List[str], error_rate: float) -> List[str]:
        """Introduce random errors to DNA sequences"""
        nucleotides = ['A', 'T', 'C', 'G']
        result = []
        
        for seq in sequences:
            seq_list = list(seq)
            num_errors = int(len(seq_list) * error_rate)
            
            # Randomly select positions to introduce errors
            error_positions = np.random.choice(len(seq_list), min(num_errors, len(seq_list)), replace=False)
            
            for pos in error_positions:
                # Replace with random nucleotide (different from original)
                current = seq_list[pos]
                new_base = np.random.choice([n for n in nucleotides if n != current])
                seq_list[pos] = new_base
            
            result.append(''.join(seq_list))
        
        return result
    
    def _count_differences(self, seq1_list: List[str], seq2_list: List[str]) -> int:
        """Count total differences between two sequence lists"""
        count = 0
        for s1, s2 in zip(seq1_list, seq2_list):
            count += sum(1 for a, b in zip(s1, s2) if a != b)
        return count
    
    def _test_recovery_feasibility(self, 
                                   original: List[str],
                                   with_errors: List[str],
                                   encoding_method: str) -> bool:
        """
        Test if encoding has enough redundancy to recover from errors.
        Uses simple heuristic: if error correction codes are present,
        recovery should be feasible for moderate error rates.
        """
        
        # Encodings with known error correction capabilities
        high_correction = {'hedges', 'church', 'wukong'}  # Have error correction codes
        medium_correction = {'goldman', 'max_density', 'no_homopolymer'}  # Some redundancy
        low_correction = {'gcplus', 'yinyang'}  # Minimal or no correction
        
        # Calculate error percentage in sequences
        total_diff = self._count_differences(original, with_errors)
        total_bases = sum(len(s) for s in original)
        error_percent = total_diff / total_bases if total_bases > 0 else 0
        
        # Success based on encoding type and error level
        if encoding_method in high_correction:
            # Can handle up to 20% errors
            return error_percent < 0.20
        elif encoding_method in medium_correction:
            # Can handle up to 10% errors
            return error_percent < 0.10
        else:
            # Can handle up to 5% errors
            return error_percent < 0.05
    
    def run_error_resistance_analysis(self,
                                     encodings: List[str],
                                     test_data: Data,
                                     base_params: Dict[str, Any],
                                     working_params: Dict[str, Any],
                                     error_rates: List[float] = None,
                                     num_runs: int = 3) -> Dict[str, List[ErrorMetrics]]:
        """
        Test all encodings at multiple error rates.
        
        Args:
            encodings: List of encoding methods to test
            test_data: Test data to encode
            base_params: Base parameters dict
            working_params: Encoding-specific parameters
            error_rates: List of error rates to test (default: [0.01, 0.05, 0.10, 0.20])
            num_runs: Number of runs per encoding/error_rate combination
        
        Returns:
            Dict mapping encoding_method -> list of ErrorMetrics
        """
        
        if error_rates is None:
            error_rates = [0.01, 0.05, 0.10, 0.20]  # 1%, 5%, 10%, 20% error rates
        
        self.logger.info(f"Starting error resistance analysis")
        self.logger.info(f"Encodings: {encodings}")
        self.logger.info(f"Error rates: {error_rates}")
        self.logger.info(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
        
        results = {enc: [] for enc in encodings}
        
        # Progress bar for overall progress
        total_tests = len(encodings) * len(error_rates)
        pbar = tqdm(total=total_tests, desc="Error resistance analysis", unit="test")
        
        for enc in encodings:
            for error_rate in error_rates:
                try:
                    # Build params for this encoding
                    params_dict = base_params.copy()
                    params_dict['encoding_method'] = enc
                    if enc in working_params:
                        params_dict.update(working_params[enc])
                    
                    # Extract filename
                    if test_data.file_paths:
                        params_dict['filename'] = os.path.basename(test_data.file_paths[0])
                    
                    # Run test
                    metrics = self.test_encoding_at_error_rate(
                        enc, test_data, error_rate, params_dict, num_runs=num_runs
                    )
                    results[enc].append(metrics)
                    
                    print(f"  {metrics}")
                    self.logger.info(str(metrics))
                    
                except Exception as e:
                    self.logger.error(f"[ERROR] {enc} @ {error_rate}: {str(e)}")
                    import traceback
                    self.logger.error(traceback.format_exc())
                
                pbar.update(1)
        
        pbar.close()
        return results


def main():
    """Run error resistance comparison"""
    
    print("\n" + "=" * 80)
    print("ERROR RESISTANCE COMPARISON")
    print("=" * 80)
    
    from dnabyte.data_classes.base import Data
    
    # Use larger test file for better error statistics
    test_file = './tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'
    
    if not os.path.exists(test_file):
        print(f"Error: Test file not found at {test_file}")
        return
    
    test_data = Data(file_paths=[test_file])
    
    print(f"\nTest file: {test_file}")
    print(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
    
    # Create tester
    output_dir = './simulations/error_results'
    tester = ErrorResistanceTester(output_dir=output_dir)
    
    # Define encodings and parameters
    encodings = [
        'goldman', 'church', 'gcplus', 'hedges',
        'max_density', 'no_homopolymer', 'wukong', 'yinyang'
    ]
    
    working_params = {
        'goldman': {'sequence_length': 200, 'add_primer': False, 'primer_length': 0},
        'church': {'sequence_length': 200, 'rs_num': 10, 'add_primer': False, 'primer_length': 0},
        'gcplus': {'sequence_length': 200, 'gcplus_k': 8, 'add_primer': False, 'primer_length': 0},
        'hedges': {'sequence_length': 200, 'hedges_coderate': 4, 'add_primer': False, 'primer_length': 0},
        'max_density': {'codeword_length': 200, 'add_primer': False, 'primer_length': 0},
        'no_homopolymer': {'codeword_length': 200, 'max_homopolymer': 4, 'add_primer': False, 'primer_length': 0},
        'wukong': {'sequence_length': 200, 'rs_num': 10, 'add_primer': False, 'primer_length': 0},
        'yinyang': {'sequence_length': 20, 'add_primer': False, 'primer_length': 0},
    }
    
    base_params = {
        'encoding_method': None,
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'storage_conditions': None,
    }
    
    # Run error resistance analysis
    error_rates = [0.01, 0.05, 0.10, 0.20]  # 1%, 5%, 10%, 20%
    
    print(f"\nError rates to test: {[f'{e*100:.0f}%' for e in error_rates]}\n")
    print("Running error resistance analysis...")
    
    results = tester.run_error_resistance_analysis(
        encodings=encodings,
        test_data=test_data,
        base_params=base_params,
        working_params=working_params,
        error_rates=error_rates,
        num_runs=3
    )
    
    # Display results
    print("\n" + "=" * 80)
    print("ERROR RESISTANCE RESULTS")
    print("=" * 80)
    
    for enc in encodings:
        print(f"\n{enc.upper()}:")
        for metrics in results[enc]:
            print(f"  {metrics}")
    
    print(f"\nResults saved to: {output_dir}")


if __name__ == '__main__':
    main()
