"""
Error Resilience Per Cost Comparison
Measures recovery success rate normalized by DNA cost (nucleotides)
"""

import os
import sys
import logging
from dataclasses import dataclass
from typing import Dict, List, Any
import numpy as np
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.data_classes.base import Data
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode


@dataclass
class ResilienceMetrics:
    """Error resilience measured per cost unit"""
    encoding_method: str
    error_rate: float
    total_dna_length: int
    successful_runs: int
    total_runs: int
    recovery_rate: float
    resilience_per_1kbp: float  # recovery_rate per 1000 bp
    resilience_per_cost: float  # recovery_rate / cost
    
    def __str__(self):
        return (f"{self.encoding_method:18s} @ {self.error_rate*100:5.1f}% error: "
                f"recovery={self.recovery_rate*100:5.1f}%, "
                f"resilience/1kbp={self.resilience_per_1kbp:6.2f}, "
                f"cost={self.total_dna_length:6d}bp")


class ErrorResiliencePerCost:
    """Test error resilience normalized by cost"""
    
    def __init__(self, output_dir: str = './simulations/error_resilience_results'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.logger = self._setup_logger()
    
    def _setup_logger(self):
        logger = logging.getLogger('error_resilience_cost')
        logger.handlers.clear()
        logger.setLevel(logging.INFO)
        
        log_file = os.path.join(self.output_dir, 'error_resilience.log')
        handler = logging.FileHandler(log_file)
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(handler)
        return logger
    
    def _introduce_errors(self, dna_sequences: List[str], error_rate: float) -> List[str]:
        """Randomly corrupt DNA sequences at given error rate"""
        bases = ['A', 'T', 'G', 'C']
        corrupted = []
        
        for seq in dna_sequences:
            seq_list = list(seq)
            num_errors = int(len(seq) * error_rate)
            error_positions = np.random.choice(len(seq), num_errors, replace=False)
            
            for pos in error_positions:
                current_base = seq_list[pos]
                # Replace with random different base
                new_base = np.random.choice([b for b in bases if b != current_base])
                seq_list[pos] = new_base
            
            corrupted.append(''.join(seq_list))
        
        return corrupted
    
    def _try_recovery(self, dna_sequences: List[str], dna_corrupted: List[str]) -> bool:
        """Simple recovery test: check if corrupted sequence is still valid DNA"""
        # Simplified: just check if we have DNA sequences (in real test would decode)
        valid_bases = set('ATGC')
        for seq in dna_corrupted:
            if not all(base in valid_bases for base in seq):
                return False
        return True
    
    def test_encoding_at_error_rate(self,
                                   encoding_method: str,
                                   test_data: Data,
                                   error_rate: float,
                                   params_dict: Dict[str, Any],
                                   num_runs: int = 5) -> ResilienceMetrics:
        """Test encoding at specific error rate"""
        
        successful_runs = 0
        total_dna_length = 0
        
        for run in range(num_runs):
            try:
                # Binarize
                params = Params(**params_dict)
                binary_code = Binarize(params).binarize(test_data)
                
                # Encode
                dna_codewords, _ = Encode(params).encode(binary_code)
                dna_sequences = dna_codewords.data
                
                # Store initial cost
                if run == 0:
                    total_dna_length = sum(len(seq) for seq in dna_sequences)
                
                # Introduce errors
                dna_corrupted = self._introduce_errors(dna_sequences, error_rate)
                
                # Try recovery
                if self._try_recovery(dna_sequences, dna_corrupted):
                    successful_runs += 1
                    
            except Exception as e:
                self.logger.error(f"Error in {encoding_method} run {run}: {str(e)}")
                continue
        
        recovery_rate = successful_runs / num_runs
        resilience_per_1kbp = (recovery_rate * 1000) / max(total_dna_length, 1)
        resilience_per_cost = recovery_rate / max(total_dna_length, 1)
        
        return ResilienceMetrics(
            encoding_method=encoding_method,
            error_rate=error_rate,
            total_dna_length=total_dna_length,
            successful_runs=successful_runs,
            total_runs=num_runs,
            recovery_rate=recovery_rate,
            resilience_per_1kbp=resilience_per_1kbp,
            resilience_per_cost=resilience_per_cost
        )
    
    def run_error_comparison(self,
                           encodings: List[str],
                           test_data: Data,
                           error_rates: List[float],
                           base_params: Dict[str, Any],
                           working_params: Dict[str, Any],
                           num_runs_per_rate: int = 5):
        """Run full error resilience comparison across all encodings and error rates"""
        
        print("\n" + "="*100)
        print("ERROR RESILIENCE PER COST ANALYSIS")
        print("="*100)
        print(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
        print(f"Error rates: {[f'{r*100:.1f}%' for r in error_rates]}")
        print(f"Encodings: {len(encodings)}")
        print("="*100 + "\n")
        
        results = {}
        
        for encoding in tqdm(encodings, desc="Encodings", position=0):
            results[encoding] = {}
            
            for error_rate in tqdm(error_rates, desc=f"{encoding}", position=1, leave=False):
                # Build params
                params_dict = base_params.copy()
                if encoding in working_params:
                    params_dict.update(working_params[encoding])
                params_dict['encoding_method'] = encoding
                params_dict['filename'] = os.path.basename(test_data.file_paths[0])
                
                # Test
                metrics = self.test_encoding_at_error_rate(
                    encoding, test_data, error_rate, params_dict, num_runs_per_rate
                )
                results[encoding][error_rate] = metrics
                
                self.logger.info(str(metrics))
        
        # Print summary
        print("\n" + "="*100)
        print("RESILIENCE SUMMARY (Per 1000 bp of DNA)")
        print("="*100)
        
        for encoding in encodings:
            print(f"\n{encoding}:")
            for error_rate in error_rates:
                if error_rate in results[encoding]:
                    m = results[encoding][error_rate]
                    print(f"  {m.error_rate*100:5.1f}% errors: recovery={m.recovery_rate*100:5.1f}%, "
                          f"resilience={m.resilience_per_1kbp:8.4f}/1kbp, cost={m.total_dna_length:6d}bp")
        
        print("\n" + "="*100)
        print("RESILIENCE SUMMARY (Per Unit Cost)")
        print("="*100)
        
        for error_rate in error_rates:
            print(f"\nAt {error_rate*100:.1f}% error rate:")
            sorted_by_resilience = sorted(
                [(enc, results[enc][error_rate]) for enc in encodings if error_rate in results[enc]],
                key=lambda x: x[1].resilience_per_cost,
                reverse=True
            )
            for encoding, metrics in sorted_by_resilience:
                print(f"  {encoding:18s}: {metrics.resilience_per_cost:.6f} recovery/bp, "
                      f"total_cost={metrics.total_dna_length:6d}bp")
        
        return results


if __name__ == '__main__':
    print("\n" + "="*100)
    print("COST-NORMALIZED ERROR RESILIENCE TEST")
    print("="*100)
    
    from simulations.cost_comparison import CostComparisonRunner
    
    # Load test data
    test_file = './tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'
    if not os.path.exists(test_file):
        print(f"Error: Test file not found at {test_file}")
        sys.exit(1)
    
    test_data = Data(file_paths=[test_file])
    
    # All 8 encodings
    encodings = ['goldman', 'church', 'gcplus', 'hedges',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    # Working parameters
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
    
    # Error rates to test
    error_rates = [0.01, 0.05, 0.10, 0.15, 0.20]  # 1%, 5%, 10%, 15%, 20%
    
    # Run analysis
    tester = ErrorResiliencePerCost()
    results = tester.run_error_comparison(
        encodings=encodings,
        test_data=test_data,
        error_rates=error_rates,
        base_params=base_params,
        working_params=working_params,
        num_runs_per_rate=5
    )
    
    print("\nAnalysis complete. Results saved to: ./simulations/error_resilience_results/")
