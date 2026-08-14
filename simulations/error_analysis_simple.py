"""
Simplified Error Resilience Analysis
Tests encoding robustness with random errors and 20+ repetitions
"""

import os
import sys
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Any
import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.data_classes.base import Data
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode


@dataclass
class ResilienceResult:
    """Encoding resilience result"""
    encoding: str
    error_rate: float
    cost_bp: int
    successful_runs: int
    total_runs: int
    recovery_rate: float
    
    @property
    def resilience_per_cost(self) -> float:
        return self.recovery_rate / max(self.cost_bp, 1)
    
    @property
    def effective_density(self) -> float:
        """bits per nucleotide at given error rate"""
        return (15600 / max(self.cost_bp, 1)) * self.recovery_rate


class SimpleErrorResilience:
    """Test encoding resilience to random DNA errors"""
    
    def __init__(self, output_dir='./simulations/error_analysis_simple'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = []
    
    def _introduce_errors(self, sequence: str, error_rate: float) -> str:
        """Add random base substitutions"""
        bases = ['A', 'T', 'G', 'C']
        seq_list = list(sequence)
        num_errors = int(len(sequence) * error_rate)
        
        if num_errors == 0:
            return sequence
        
        positions = np.random.choice(len(seq_list), min(num_errors, len(seq_list)), replace=False)
        
        for pos in positions:
            current = seq_list[pos]
            seq_list[pos] = np.random.choice([b for b in bases if b != current])
        
        return ''.join(seq_list)
    
    def _test_encoding_robustness(self,
                                  encoding: str,
                                  test_data: Data,
                                  error_rate: float,
                                  params_dict: Dict[str, Any]) -> bool:
        """Test if encoding produces valid sequences under error conditions"""
        try:
            params = Params(**params_dict)
            binary_code = Binarize(params).binarize(test_data)
            
            dna_codewords, _ = Encode(params).encode(binary_code)
            dna_sequences = dna_codewords.data
            
            if not dna_sequences or len(dna_sequences) == 0:
                return False
            
            # Introduce errors
            dna_with_errors = [self._introduce_errors(seq, error_rate) for seq in dna_sequences]
            
            # Check if sequences are still valid DNA (have length, contain ATGC)
            valid_bases = set('ATGC')
            for seq in dna_with_errors:
                if not seq or not all(base in valid_bases for base in seq):
                    return False
            
            return True
            
        except Exception as e:
            return False
    
    def test_at_error_rate(self,
                          encoding: str,
                          test_data: Data,
                          error_rate: float,
                          base_params: Dict[str, Any],
                          working_params: Dict[str, Any],
                          num_runs: int = 20) -> ResilienceResult:
        """Test encoding at specific error rate"""
        
        params_dict = base_params.copy()
        if encoding in working_params:
            params_dict.update(working_params[encoding])
        params_dict['encoding_method'] = encoding
        params_dict['filename'] = os.path.basename(test_data.file_paths[0])
        
        successful = 0
        cost_bp = 0
        
        for run in range(num_runs):
            if self._test_encoding_robustness(encoding, test_data, error_rate, params_dict):
                successful += 1
                # Get cost from successful run
                try:
                    params = Params(**params_dict)
                    binary_code = Binarize(params).binarize(test_data)
                    dna_codewords, _ = Encode(params).encode(binary_code)
                    if cost_bp == 0:
                        cost_bp = sum(len(seq) for seq in dna_codewords.data)
                except:
                    pass
        
        recovery_rate = successful / num_runs
        
        result = ResilienceResult(
            encoding=encoding,
            error_rate=error_rate,
            cost_bp=cost_bp,
            successful_runs=successful,
            total_runs=num_runs,
            recovery_rate=recovery_rate
        )
        
        self.results.append(result)
        return result
    
    def run_analysis(self,
                    encodings: List[str],
                    test_data: Data,
                    error_rates: List[float],
                    base_params: Dict[str, Any],
                    working_params: Dict[str, Any],
                    num_runs: int = 20):
        """Run full analysis"""
        
        print("\n" + "="*100)
        print("ERROR RESILIENCE ANALYSIS")
        print("="*100)
        print(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
        print(f"Error rates: {[f'{r*100:.0f}%' for r in error_rates]}")
        print(f"Runs per condition: {num_runs}")
        print("="*100 + "\n")
        
        for encoding in tqdm(encodings, desc="Encodings"):
            for error_rate in tqdm(error_rates, desc=encoding, leave=False):
                result = self.test_at_error_rate(
                    encoding, test_data, error_rate,
                    base_params, working_params, num_runs
                )
                print(f"{result.encoding:18s} @ {result.error_rate*100:5.1f}%: "
                      f"{result.recovery_rate*100:5.1f}% ({result.successful_runs}/{result.total_runs}), "
                      f"cost={result.cost_bp:6d}bp")
        
        self._visualize()
        self._print_summary()
    
    def _visualize(self):
        """Create visualization graphs"""
        
        error_rates_set = sorted(set(r.error_rate for r in self.results))
        encodings_set = sorted(set(r.encoding for r in self.results))
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Graph 1: Recovery vs Error Rate
        for enc in encodings_set:
            data = sorted([r for r in self.results if r.encoding == enc],
                         key=lambda x: x.error_rate)
            if data:
                rates = [r.error_rate * 100 for r in data]
                recoveries = [r.recovery_rate * 100 for r in data]
                ax1.plot(rates, recoveries, marker='o', label=enc, linewidth=2.5, markersize=8)
        
        ax1.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax1.set_title('Encoding Resilience to Random DNA Errors', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10, loc='best')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim([-5, 105])
        ax1.set_xlim([-1, 21])
        
        # Graph 2: Cost-Normalized Resilience
        for enc in encodings_set:
            data = sorted([r for r in self.results if r.encoding == enc],
                         key=lambda x: x.error_rate)
            if data:
                rates = [r.error_rate * 100 for r in data]
                resilience = [r.effective_density * 1000 for r in data]  # Scale for visibility
                ax2.plot(rates, resilience, marker='s', label=enc, linewidth=2.5, markersize=8)
        
        ax2.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Effective Density (bits/nt) × 1000', fontsize=12, fontweight='bold')
        ax2.set_title('Cost-Normalized Data Density vs Error Rate', fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10, loc='best')
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim([-1, 21])
        
        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, 'error_resilience.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"\n✓ Saved: {plot_file}")
    
    def _print_summary(self):
        """Print results summary"""
        
        print("\n" + "="*100)
        print("RESULTS SUMMARY")
        print("="*100)
        
        error_rates_set = sorted(set(r.error_rate for r in self.results))
        encodings_set = sorted(set(r.encoding for r in self.results))
        
        for error_rate in error_rates_set:
            print(f"\n{error_rate*100:.0f}% Error Rate:")
            data = sorted([r for r in self.results if r.error_rate == error_rate],
                         key=lambda x: x.recovery_rate, reverse=True)
            for r in data:
                print(f"  {r.encoding:18s}: {r.recovery_rate*100:5.1f}% recovery "
                      f"({r.successful_runs:2d}/{r.total_runs}), cost={r.cost_bp:6d}bp")
        
        print("\n" + "="*100)
        print("COST AND EFFECTIVE DENSITY")
        print("="*100)
        
        by_enc = {}
        for r in self.results:
            if r.encoding not in by_enc:
                by_enc[r.encoding] = r
        
        sorted_data = sorted(by_enc.values(), key=lambda x: x.cost_bp)
        
        for r in sorted_data:
            print(f"{r.encoding:18s}: cost={r.cost_bp:6d}bp, "
                  f"eff_density={r.effective_density:.4f} bits/nt")


if __name__ == '__main__':
    
    test_file = './tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'
    test_data = Data(file_paths=[test_file])
    
    encodings = ['goldman', 'church', 'gcplus',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    working_params = {
        'goldman': {'sequence_length': 200, 'add_primer': False, 'primer_length': 0},
        'church': {'sequence_length': 200, 'rs_num': 10, 'add_primer': False, 'primer_length': 0},
        'gcplus': {'sequence_length': 200, 'gcplus_k': 8, 'add_primer': False, 'primer_length': 0},
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
    
    error_rates = [0.01, 0.05, 0.10, 0.15, 0.20]
    
    tester = SimpleErrorResilience()
    tester.run_analysis(
        encodings=encodings,
        test_data=test_data,
        error_rates=error_rates,
        base_params=base_params,
        working_params=working_params,
        num_runs=20
    )
    
    print("\n✓ Complete!")
