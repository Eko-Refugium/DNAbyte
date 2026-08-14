"""
Real Error Resilience Analysis
Tests encoding recovery using actual Cluster + Consensus pipeline with injected errors
"""

import os
import sys
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Any, Tuple
import matplotlib.pyplot as plt
from tqdm import tqdm
import logging

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.data_classes.base import Data
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus


@dataclass
class ResilienceResult:
    """Encoding resilience result"""
    encoding: str
    error_rate: float
    cost_bp: int
    successful_runs: int
    total_runs: int
    
    @property
    def recovery_rate(self) -> float:
        return self.successful_runs / self.total_runs
    
    @property
    def resilience_per_cost(self) -> float:
        return self.recovery_rate / max(self.cost_bp, 1)
    
    @property
    def effective_density(self) -> float:
        """bits per nucleotide at given error rate"""
        return (15600 / max(self.cost_bp, 1)) * self.recovery_rate


class RealErrorResilience:
    """Test encoding recovery with actual Cluster+Consensus pipeline"""
    
    def __init__(self, output_dir='./simulations/error_analysis_real'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = []
        
        # Set up logging
        log_file = os.path.join(output_dir, 'recovery_test.log')
        logging.basicConfig(
            filename=log_file,
            level=logging.INFO,
            format='%(message)s',
            filemode='w'
        )
    
    def _introduce_errors(self, sequence: str, error_rate: float) -> str:
        """Add random base substitutions to sequence"""
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
    
    def _test_single_recovery(self,
                             encoding: str,
                             binary_code,
                             params_dict: Dict[str, Any],
                             error_rate: float) -> Tuple[bool, int]:
        """
        Test recovery of data through actual Cluster+Consensus pipeline
        Returns: (success, cost_bp)
        """
        try:
            # Encode original data
            params = Params(**params_dict)
            dna_codewords, _ = Encode(params).encode(binary_code)
            dna_sequences = dna_codewords.data
            
            if not dna_sequences or len(dna_sequences) == 0:
                return False, 0
            
            cost_bp = sum(len(seq) for seq in dna_sequences)
            
            # Introduce errors
            dna_with_errors = [self._introduce_errors(seq, error_rate) for seq in dna_sequences]
            
            # Create InSilicoDNA object with errored sequences
            dna_obj = InSilicoDNA(dna_with_errors)
            
            # Try to recover using Cluster + Consensus
            cluster_obj = Cluster(params)
            try:
                clustered_data, _ = cluster_obj.cluster(dna_obj)
            except Exception as e:
                return False, cost_bp
            
            consensus_obj = Consensus(params)
            try:
                recovered_data, _ = consensus_obj.call(clustered_data)
            except Exception as e:
                return False, cost_bp
            
            # Success = recovered data has content
            success = len(recovered_data.data) > 0
            return success, cost_bp
            
        except Exception as e:
            return False, 0
    
    def test_encoding_at_error_rate(self,
                                   encoding: str,
                                   test_data: Data,
                                   error_rate: float,
                                   base_params: Dict[str, Any],
                                   working_params: Dict[str, Any],
                                   num_runs: int = 20) -> ResilienceResult:
        """Test encoding at specific error rate with real recovery pipeline"""
        
        # Binarize once (same for all runs)
        params_dict = base_params.copy()
        if encoding in working_params:
            params_dict.update(working_params[encoding])
        params_dict['encoding_method'] = encoding
        params_dict['filename'] = os.path.basename(test_data.file_paths[0])
        
        params = Params(**params_dict)
        binary_code = Binarize(params).binarize(test_data)
        
        successful = 0
        cost_bp = 0
        
        for run in range(num_runs):
            success, bp = self._test_single_recovery(encoding, binary_code, params_dict, error_rate)
            if success:
                successful += 1
            if cost_bp == 0:
                cost_bp = bp
        
        recovery_rate = successful / num_runs
        
        result = ResilienceResult(
            encoding=encoding,
            error_rate=error_rate,
            cost_bp=cost_bp,
            successful_runs=successful,
            total_runs=num_runs
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
        """Run full analysis with real recovery pipeline"""
        
        print("\n" + "="*100)
        print("REAL ERROR RESILIENCE ANALYSIS")
        print("="*100)
        print(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
        print(f"Error rates: {[f'{r*100:.0f}%' for r in error_rates]}")
        print(f"Runs per condition: {num_runs}")
        print("Recovery method: Cluster + Consensus (actual pipeline)")
        print("="*100 + "\n")
        
        for encoding in tqdm(encodings, desc="Encodings"):
            for error_rate in tqdm(error_rates, desc=encoding, leave=False):
                try:
                    result = self.test_encoding_at_error_rate(
                        encoding, test_data, error_rate,
                        base_params, working_params, num_runs
                    )
                    print(f"{result.encoding:18s} @ {result.error_rate*100:5.1f}%: "
                          f"{result.recovery_rate*100:5.1f}% ({result.successful_runs}/{result.total_runs}), "
                          f"cost={result.cost_bp:6d}bp")
                except Exception as e:
                    print(f"{encoding:18s} @ {error_rate*100:5.1f}%: ERROR - {str(e)[:50]}")
        
        self._visualize()
        self._print_summary()
    
    def _visualize(self):
        """Create visualization graphs"""
        
        if not self.results:
            print("No results to visualize")
            return
        
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
        ax1.set_title('Real Recovery: Cluster+Consensus vs Error Rate', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10, loc='best')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim([-5, 105])
        ax1.set_xlim([-1, 21])
        
        # Graph 2: Cost vs Recovery (scatter)
        for error_rate in error_rates_set:
            data = [r for r in self.results if r.error_rate == error_rate]
            if data:
                costs = [r.cost_bp / 1000 for r in data]
                recoveries = [r.recovery_rate * 100 for r in data]
                ax2.scatter(costs, recoveries, s=150, alpha=0.7, 
                           label=f'{error_rate*100:.0f}% error')
        
        ax2.set_xlabel('Cost (thousands of bp)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax2.set_title('Recovery vs Cost (colored by error rate)', fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10, loc='best')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, 'recovery_analysis.png')
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
        print("BEST ENCODING BY ERROR RATE")
        print("="*100)
        
        for error_rate in error_rates_set:
            data = sorted([r for r in self.results if r.error_rate == error_rate],
                         key=lambda x: x.recovery_rate, reverse=True)
            if data:
                best = data[0]
                print(f"{error_rate*100:5.1f}%: {best.encoding:18s} ({best.recovery_rate*100:5.1f}%)")


if __name__ == '__main__':
    
    test_file = './tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'
    test_data = Data(file_paths=[test_file])
    
    encodings = ['goldman', 'church', 'gcplus',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    working_params = {
        'goldman': {
            'sequence_length': 200,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'church': {
            'sequence_length': 200,
            'rs_num': 10,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'gcplus': {
            'sequence_length': 200,
            'gcplus_k': 8,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'max_density': {
            'codeword_length': 200,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'no_homopolymer': {
            'codeword_length': 200,
            'max_homopolymer': 4,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'wukong': {
            'sequence_length': 200,
            'rs_num': 10,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'yinyang': {
            'sequence_length': 20,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
    }
    
    base_params = {
        'encoding_method': None,
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'storage_conditions': None,
    }
    
    error_rates = [0.01, 0.05, 0.10, 0.15, 0.20]
    
    tester = RealErrorResilience()
    tester.run_analysis(
        encodings=encodings,
        test_data=test_data,
        error_rates=error_rates,
        base_params=base_params,
        working_params=working_params,
        num_runs=20
    )
    
    print("\n✓ Complete!")
