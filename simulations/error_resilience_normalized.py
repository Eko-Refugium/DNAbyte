"""
Cost-Normalized Error Resilience Analysis
All encodings with sequence_length=100bp, cost normalized via mean parameter
Uses full pipeline with kmere sequencing errors
"""

import os
import sys
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Any
import matplotlib.pyplot as plt
from tqdm import tqdm
import logging

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.data_classes.base import Data
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.sequence import SimulateSequencing
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus
from dnabyte.data_classes.insilicodna import InSilicoDNA


@dataclass
class ResilienceResult:
    """Encoding resilience result"""
    encoding: str
    error_rate: float
    cost_bp: int
    mean: int
    successful_runs: int
    total_runs: int
    
    @property
    def recovery_rate(self) -> float:
        return self.successful_runs / self.total_runs
    
    @property
    def resilience_per_cost(self) -> float:
        return self.recovery_rate / max(self.cost_bp, 1)


class CostNormalizedResilience:
    """Test encodings with cost-normalized parameters"""
    
    def __init__(self, output_dir='./simulations/error_resilience_normalized'):
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
    
    def _test_single_run(self,
                        encoding: str,
                        binary_code,
                        params_dict: Dict[str, Any],
                        error_rate: float) -> tuple:
        """
        Full pipeline test: Encode -> Synthesize -> Introduce errors -> Sequence -> Cluster -> Consensus
        Returns: (success, cost_bp)
        """
        try:
            params = Params(**params_dict)
            
            # Encode
            enc_obj = Encode(params)
            dna_codewords, _ = enc_obj.encode(binary_code)
            dna_sequences = dna_codewords.data
            
            if not dna_sequences or len(dna_sequences) == 0:
                return False, 0
            
            # Get baseline cost
            baseline_cost = sum(len(seq) for seq in dna_sequences)
            
            # Synthesize (create copies)
            syn_obj = SimulateSynthesis(params)
            dna_synthesized, _ = syn_obj.simulate(dna_codewords)
            dna_syn_list = dna_synthesized.data
            
            # Introduce random errors
            dna_with_errors = [self._introduce_errors(seq, error_rate) for seq in dna_syn_list]
            
            # Simulate sequencing
            dna_obj = InSilicoDNA(dna_with_errors)
            
            seq_obj = SimulateSequencing(params)
            dna_sequenced, _ = seq_obj.simulate(dna_obj)
            
            # Cluster and consensus
            cluster_obj = Cluster(params)
            clustered_data, _ = cluster_obj.cluster(dna_sequenced)
            
            consensus_obj = Consensus(params)
            recovered_data, _ = consensus_obj.call(clustered_data)
            
            # Success if we recovered data
            success = len(recovered_data.data) > 0
            return success, baseline_cost
            
        except Exception as e:
            return False, 0
    
    def calculate_mean_for_target_cost(self,
                                      encoding: str,
                                      test_data: Data,
                                      target_cost_bp: int,
                                      base_params: Dict[str, Any],
                                      encoding_params: Dict[str, Any]) -> int:
        """Calculate mean value needed to achieve target cost"""
        
        # Test with mean=1 to get baseline
        params_dict = base_params.copy()
        params_dict.update(encoding_params)
        params_dict['encoding_method'] = encoding
        params_dict['mean'] = 1
        params_dict['filename'] = os.path.basename(test_data.file_paths[0])
        
        try:
            params = Params(**params_dict)
            binary_code = Binarize(params).binarize(test_data)
            dna_codewords, _ = Encode(params).encode(binary_code)
            baseline_cost = sum(len(seq) for seq in dna_codewords.data)
            
            # Calculate mean needed
            mean_needed = target_cost_bp / baseline_cost
            return max(1, int(mean_needed))
        except:
            return 1
    
    def run_analysis(self,
                    encodings: List[str],
                    test_data: Data,
                    target_cost_bp: int,
                    error_rates: List[float],
                    base_params: Dict[str, Any],
                    encoding_params: Dict[str, Dict[str, Any]],
                    num_runs: int = 5):
        """Run cost-normalized analysis"""
        
        print("\n" + "="*100)
        print("COST-NORMALIZED ERROR RESILIENCE ANALYSIS")
        print("="*100)
        print(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
        print(f"Sequence length: 100 bp (fixed)")
        print(f"Target cost: {target_cost_bp:,} bp (equal for all)")
        print(f"Error rates: {[f'{r*100:.0f}%' for r in error_rates]}")
        print(f"Runs per condition: {num_runs}")
        print(f"Recovery: Full pipeline (Encode->Synthesize->Sequence->Cluster->Consensus)")
        print("="*100 + "\n")
        
        # Step 1: Calculate mean values for each encoding
        print("Step 1: Calculating mean values to normalize costs...\n")
        mean_values = {}
        
        for encoding in encodings:
            mean_val = self.calculate_mean_for_target_cost(
                encoding, test_data, target_cost_bp, base_params, encoding_params[encoding]
            )
            mean_values[encoding] = mean_val
            print(f"  {encoding:18s}: mean={mean_val:3d} (target cost=[OK] {target_cost_bp:,} bp)")
        
        # Step 2: Run error resilience tests
        print("\n" + "="*100)
        print("Step 2: Running error resilience tests...")
        print("="*100 + "\n")
        
        for encoding in tqdm(encodings, desc="Encodings"):
            binary_code = None
            
            for error_rate in tqdm(error_rates, desc=encoding, leave=False):
                # Binarize once per encoding
                if binary_code is None:
                    params_dict = base_params.copy()
                    params_dict.update(encoding_params[encoding])
                    params_dict['encoding_method'] = encoding
                    params_dict['filename'] = os.path.basename(test_data.file_paths[0])
                    params = Params(**params_dict)
                    binary_code = Binarize(params).binarize(test_data)
                
                # Test with normalized cost
                params_dict = base_params.copy()
                params_dict.update(encoding_params[encoding])
                params_dict['encoding_method'] = encoding
                params_dict['mean'] = mean_values[encoding]
                params_dict['filename'] = os.path.basename(test_data.file_paths[0])
                
                successful = 0
                cost_bp = 0
                
                for run in range(num_runs):
                    success, bp = self._test_single_run(encoding, binary_code, params_dict, error_rate)
                    if success:
                        successful += 1
                    if cost_bp == 0:
                        cost_bp = bp * mean_values[encoding]  # Total cost = baseline * mean
                
                result = ResilienceResult(
                    encoding=encoding,
                    error_rate=error_rate,
                    cost_bp=cost_bp,
                    mean=mean_values[encoding],
                    successful_runs=successful,
                    total_runs=num_runs
                )
                self.results.append(result)
                
                print(f"{encoding:18s} @ {error_rate*100:5.1f}%: "
                      f"{result.recovery_rate*100:5.1f}% ({successful}/{num_runs}), "
                      f"cost={cost_bp:7d}bp, mean={mean_values[encoding]:3d}")
        
        self._visualize()
        self._print_summary()
    
    def _visualize(self):
        """Create visualization"""
        
        if not self.results:
            return
        
        error_rates_set = sorted(set(r.error_rate for r in self.results))
        encodings_set = sorted(set(r.encoding for r in self.results))
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Recovery vs Error Rate (all costs equal)
        for enc in encodings_set:
            data = sorted([r for r in self.results if r.encoding == enc],
                         key=lambda x: x.error_rate)
            if data:
                rates = [r.error_rate * 100 for r in data]
                recoveries = [r.recovery_rate * 100 for r in data]
                ax.plot(rates, recoveries, marker='o', label=enc, linewidth=2.5, markersize=8)
        
        ax.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Recovery vs Error Rate (All Encodings at Equal Cost)', 
                     fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([-5, 105])
        ax.set_xlim([-1, 21])
        
        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, 'cost_normalized_recovery.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"\n[OK] Saved: {plot_file}")
    
    def _print_summary(self):
        """Print summary"""
        
        print("\n" + "="*100)
        print("SUMMARY")
        print("="*100)
        
        error_rates_set = sorted(set(r.error_rate for r in self.results))
        encodings_set = sorted(set(r.encoding for r in self.results))
        
        for error_rate in error_rates_set:
            print(f"\n{error_rate*100:.0f}% Error Rate:")
            data = sorted([r for r in self.results if r.error_rate == error_rate],
                         key=lambda x: x.recovery_rate, reverse=True)
            for r in data:
                print(f"  {r.encoding:18s}: {r.recovery_rate*100:5.1f}% recovery "
                      f"({r.successful_runs:2d}/{r.total_runs}), mean={r.mean:3d}")


if __name__ == '__main__':
    
    test_file = './tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'
    test_data = Data(file_paths=[test_file])
    
    encodings = ['goldman', 'church', 'gcplus', 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    # All with sequence_length=100bp (fixed)
    encoding_params = {
        'goldman': {'sequence_length': 100, 'add_primer': False, 'primer_length': 0},
        'church': {'sequence_length': 100, 'rs_num': 10, 'add_primer': False, 'primer_length': 0},
        'gcplus': {'sequence_length': 100, 'gcplus_k': 8, 'add_primer': False, 'primer_length': 0},
        'max_density': {'codeword_length': 100, 'add_primer': False, 'primer_length': 0},
        'no_homopolymer': {'codeword_length': 100, 'max_homopolymer': 4, 'add_primer': False, 'primer_length': 0},
        'wukong': {'sequence_length': 100, 'rs_num': 10, 'add_primer': False, 'primer_length': 0},
        'yinyang': {'sequence_length': 100, 'add_primer': False, 'primer_length': 0},
    }
    
    base_params = {
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'kmere',
        'kmer_k': 1,
        'kmer_p_ins': 0.001,
        'kmer_p_del': 0.001,
        'kmer_p_sub': 0.001,
        'kmer_seed': 42,
        'storage_conditions': None,
    }
    
    error_rates = [0.01, 0.05, 0.10, 0.15, 0.20]
    target_cost = 20000  # bp - all encodings normalized to this
    
    tester = CostNormalizedResilience()
    tester.run_analysis(
        encodings=encodings,
        test_data=test_data,
        target_cost_bp=target_cost,
        error_rates=error_rates,
        base_params=base_params,
        encoding_params=encoding_params,
        num_runs=5
    )
    
    print("\n[OK] Complete!")
    )
    
    print("\n[OK] Complete!")
