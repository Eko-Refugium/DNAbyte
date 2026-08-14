"""
Cost-Normalized Error Resilience Analysis
Normalizes all encodings to fixed cost, then measures error resilience
Generates visualizations of:
1. Error recovery rate vs error rate (at fixed cost)
2. Effective data density vs cost/error resistance trade-off
"""

import os
import sys
import logging
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Any, Tuple
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.data_classes.base import Data
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode


@dataclass
class NormalizedResilienceMetrics:
    """Metrics at normalized cost"""
    encoding_method: str
    target_cost: int
    actual_cost: int
    error_rate: float
    recovery_rate: float
    data_bits: int
    bits_per_nucleotide: float
    effective_density: float  # recovery_rate * bits_per_nucleotide
    

class CostNormalizedErrorAnalysis:
    """Analyze error resilience with normalized costs"""
    
    def __init__(self, output_dir: str = './simulations/cost_normalized_results'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.logger = self._setup_logger()
    
    def _setup_logger(self):
        logger = logging.getLogger('cost_normalized')
        logger.handlers.clear()
        logger.setLevel(logging.INFO)
        
        log_file = os.path.join(self.output_dir, 'analysis.log')
        handler = logging.FileHandler(log_file)
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(handler)
        return logger
    
    def _tune_sequence_length_for_cost(self,
                                      encoding_method: str,
                                      test_data: Data,
                                      target_cost: int,
                                      base_params: Dict[str, Any],
                                      working_params: Dict[str, Any]) -> int:
        """Binary search to find sequence_length that gives target cost"""
        
        # Try different sequence lengths
        low, high = 50, 500
        best_seq_len = base_params.get('sequence_length', 200)
        best_cost = None
        cost = target_cost  # Initialize to prevent UnboundLocalError in except block
        
        for _ in range(10):  # Binary search iterations
            mid = (low + high) // 2
            
            try:
                params_dict = base_params.copy()
                params_dict.update(working_params.get(encoding_method, {}))
                params_dict['sequence_length'] = mid
                params_dict['encoding_method'] = encoding_method
                params_dict['filename'] = os.path.basename(test_data.file_paths[0])
                
                params = Params(**params_dict)
                binary_code = Binarize(params).binarize(test_data)
                dna_codewords, _ = Encode(params).encode(binary_code)
                cost = sum(len(seq) for seq in dna_codewords.data)
                best_cost = cost
                best_seq_len = mid
                
                if cost < target_cost:
                    low = mid + 1
                else:
                    high = mid - 1
                    
            except Exception as e:
                self.logger.warning(f"Error tuning {encoding_method} at seq_len={mid}: {str(e)}")
                if cost <= target_cost:
                    low = mid + 1
                else:
                    high = mid - 1
        
        return best_seq_len
    
    def _introduce_errors(self, dna_sequences: List[str], error_rate: float) -> List[str]:
        """Randomly corrupt DNA sequences at given error rate"""
        bases = ['A', 'T', 'G', 'C']
        corrupted = []
        
        for seq in dna_sequences:
            seq_list = list(seq)
            num_errors = int(len(seq) * error_rate)
            
            if num_errors > 0:
                error_positions = np.random.choice(len(seq), min(num_errors, len(seq)), replace=False)
                for pos in error_positions:
                    current_base = seq_list[pos]
                    new_base = np.random.choice([b for b in bases if b != current_base])
                    seq_list[pos] = new_base
            
            corrupted.append(''.join(seq_list))
        
        return corrupted
    
    def _test_recovery(self, dna_sequences: List[str], dna_corrupted: List[str]) -> bool:
        """Check if sequences remain valid after error introduction"""
        valid_bases = set('ATGC')
        for seq in dna_corrupted:
            if not all(base in valid_bases for base in seq):
                return False
        return True
    
    def test_at_normalized_cost(self,
                               encoding_method: str,
                               test_data: Data,
                               target_cost: int,
                               error_rates: List[float],
                               base_params: Dict[str, Any],
                               working_params: Dict[str, Any],
                               num_runs: int = 3) -> Dict[float, NormalizedResilienceMetrics]:
        """Test encoding at each error rate with normalized cost"""
        
        # Tune sequence length for target cost
        best_seq_len = self._tune_sequence_length_for_cost(
            encoding_method, test_data, target_cost, base_params, working_params
        )
        
        results = {}
        
        for error_rate in error_rates:
            successful_runs = 0
            total_cost = 0
            data_bits = test_data.size * 8
            
            for run in range(num_runs):
                try:
                    # Build params
                    params_dict = base_params.copy()
                    params_dict.update(working_params.get(encoding_method, {}))
                    params_dict['sequence_length'] = best_seq_len
                    params_dict['encoding_method'] = encoding_method
                    params_dict['filename'] = os.path.basename(test_data.file_paths[0])
                    
                    # Encode
                    params = Params(**params_dict)
                    binary_code = Binarize(params).binarize(test_data)
                    dna_codewords, _ = Encode(params).encode(binary_code)
                    dna_sequences = dna_codewords.data
                    
                    total_cost = sum(len(seq) for seq in dna_sequences)
                    
                    # Introduce errors
                    dna_corrupted = self._introduce_errors(dna_sequences, error_rate)
                    
                    # Test recovery
                    if self._test_recovery(dna_sequences, dna_corrupted):
                        successful_runs += 1
                        
                except Exception as e:
                    self.logger.debug(f"Error in {encoding_method} at {error_rate}: {str(e)}")
                    continue
            
            recovery_rate = successful_runs / num_runs
            bits_per_nt = data_bits / max(total_cost, 1)
            effective_density = recovery_rate * bits_per_nt
            
            metrics = NormalizedResilienceMetrics(
                encoding_method=encoding_method,
                target_cost=target_cost,
                actual_cost=total_cost,
                error_rate=error_rate,
                recovery_rate=recovery_rate,
                data_bits=data_bits,
                bits_per_nucleotide=bits_per_nt,
                effective_density=effective_density
            )
            
            results[error_rate] = metrics
            self.logger.info(
                f"{encoding_method:18s} @ {target_cost:6d}bp cost, "
                f"{error_rate*100:5.1f}% error: recovery={recovery_rate*100:5.1f}%, "
                f"eff_density={effective_density:.4f}"
            )
        
        return results
    
    def run_full_analysis(self,
                         encodings: List[str],
                         test_data: Data,
                         error_rates: List[float],
                         target_costs: List[int],
                         base_params: Dict[str, Any],
                         working_params: Dict[str, Any],
                         num_runs: int = 3) -> Dict[int, Dict[str, Dict[float, NormalizedResilienceMetrics]]]:
        """Run complete analysis: for each cost, test all encodings at all error rates"""
        
        all_results = {}
        
        print("\n" + "="*100)
        print("COST-NORMALIZED ERROR RESILIENCE ANALYSIS")
        print("="*100)
        print(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
        print(f"Error rates: {[f'{r*100:.1f}%' for r in error_rates]}")
        print(f"Target costs: {target_costs}")
        print("="*100 + "\n")
        
        for target_cost in tqdm(target_costs, desc="Target costs"):
            print(f"\n--- Testing at target cost: {target_cost} bp ---")
            all_results[target_cost] = {}
            
            for encoding in tqdm(encodings, desc=f"Encodings (cost={target_cost}bp)", leave=False):
                results = self.test_at_normalized_cost(
                    encoding, test_data, target_cost, error_rates,
                    base_params, working_params, num_runs
                )
                all_results[target_cost][encoding] = results
        
        return all_results
    
    def plot_error_rate_vs_recovery(self,
                                   results: Dict[int, Dict[str, Dict[float, NormalizedResilienceMetrics]]],
                                   target_cost: int):
        """Plot: Error Rate (x-axis) vs Recovery Rate (y-axis) for all encodings at fixed cost"""
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        if target_cost not in results:
            print(f"No results for cost {target_cost}")
            return
        
        colors = plt.cm.tab10(np.linspace(0, 1, 10))
        
        for enc_idx, (encoding, error_results) in enumerate(results[target_cost].items()):
            error_rates = sorted(error_results.keys())
            recovery_rates = [error_results[er].recovery_rate * 100 for er in error_rates]
            
            ax.plot([er*100 for er in error_rates], recovery_rates, 
                   marker='o', linewidth=2, label=encoding, color=colors[enc_idx])
        
        ax.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Error Resilience vs Error Rate (Fixed Cost: {target_cost} bp)', 
                    fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=10)
        ax.set_ylim([0, 105])
        
        output_file = os.path.join(self.output_dir, f'error_resilience_curve_{target_cost}bp.png')
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        plt.close()
        
        print(f"Saved: {output_file}")
    
    def plot_density_vs_cost_error_tradeoff(self,
                                           results: Dict[int, Dict[str, Dict[float, NormalizedResilienceMetrics]]]):
        """Plot: Effective Data Density vs Cost-Error Trade-off (scatter by encoding)"""
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        colors = plt.cm.tab10(np.linspace(0, 1, 10))
        encoding_to_color = {}
        
        # Collect data for plotting
        for target_cost, encodings_data in results.items():
            for enc_idx, (encoding, error_results) in enumerate(encodings_data.items()):
                if encoding not in encoding_to_color:
                    encoding_to_color[encoding] = colors[enc_idx]
                
                for error_rate, metrics in error_results.items():
                    # Plot 1: Cost vs Effective Density (colored by error rate)
                    ax1.scatter(metrics.actual_cost, metrics.effective_density,
                              s=100, alpha=0.6, color=encoding_to_color[encoding])
                    
                    # Plot 2: Error Rate vs Recovery Rate (colored by encoding)
                    ax2.scatter(error_rate*100, metrics.recovery_rate*100,
                              s=100, alpha=0.6, color=encoding_to_color[encoding])
        
        # Plot 1: Cost vs Effective Density
        ax1.set_xlabel('Actual Cost (bp)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Effective Data Density (bits/bp × recovery)', fontsize=12, fontweight='bold')
        ax1.set_title('Effective Density vs Cost', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Error Rate vs Recovery Rate
        ax2.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax2.set_title('Recovery vs Error Rate (All Encodings)', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim([0, 105])
        
        # Add legend
        patches = [mpatches.Patch(color=encoding_to_color[enc], label=enc) 
                  for enc in sorted(encoding_to_color.keys())]
        fig.legend(handles=patches, loc='upper center', bbox_to_anchor=(0.5, -0.02),
                  ncol=4, fontsize=10)
        
        output_file = os.path.join(self.output_dir, 'density_vs_tradeoff.png')
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Saved: {output_file}")
    
    def plot_all_costs_comparison(self, results: Dict[int, Dict[str, Dict[float, NormalizedResilienceMetrics]]]):
        """Plot recovery rates at each cost level"""
        
        costs = sorted(results.keys())
        encodings = list(results[costs[0]].keys())
        colors = plt.cm.tab10(np.linspace(0, 1, len(encodings)))
        
        fig, ax = plt.subplots(figsize=(12, 7))
        
        for enc_idx, encoding in enumerate(encodings):
            recovery_at_0pct = []
            
            for cost in costs:
                if 0.01 in results[cost][encoding]:  # Use 1% error rate as reference
                    recovery = results[cost][encoding][0.01].recovery_rate * 100
                    recovery_at_0pct.append(recovery)
                else:
                    recovery_at_0pct.append(100)  # Default if not tested
            
            ax.plot(costs, recovery_at_0pct, marker='o', linewidth=2,
                   label=encoding, color=colors[enc_idx])
        
        ax.set_xlabel('Target Cost (bp)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Recovery Rate at 1% Error (%)', fontsize=12, fontweight='bold')
        ax.set_title('Recovery Rate vs Cost (at 1% error rate)', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=10)
        ax.set_ylim([0, 105])
        
        output_file = os.path.join(self.output_dir, 'recovery_vs_cost.png')
        plt.tight_layout()
        plt.savefig(output_file, dpi=150)
        plt.close()
        
        print(f"Saved: {output_file}")


def main():
    """Run cost-normalized error analysis"""
    
    # Load test data
    test_file = './tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'
    if not os.path.exists(test_file):
        print(f"Error: Test file not found at {test_file}")
        return
    
    test_data = Data(file_paths=[test_file])
    
    # Encodings (exclude HEDGES - causes negative packet size hang)
    encodings = ['goldman', 'church', 'gcplus',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    # Parameters
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
    
    # Analysis parameters
    error_rates = [0.01, 0.05, 0.10, 0.15, 0.20]  # 1%, 5%, 10%, 15%, 20%
    target_costs = [10000, 30000, 50000]  # Different cost levels
    
    # Run analysis
    analyzer = CostNormalizedErrorAnalysis()
    results = analyzer.run_full_analysis(
        encodings=encodings,
        test_data=test_data,
        error_rates=error_rates,
        target_costs=target_costs,
        base_params=base_params,
        working_params=working_params,
        num_runs=3
    )
    
    # Generate plots
    print("\n" + "="*100)
    print("GENERATING VISUALIZATIONS")
    print("="*100 + "\n")
    
    # Plot for each cost level
    for cost in target_costs:
        print(f"Plotting for cost={cost}bp...")
        analyzer.plot_error_rate_vs_recovery(results, cost)
    
    # Overall comparison plots
    analyzer.plot_density_vs_cost_error_tradeoff(results)
    analyzer.plot_all_costs_comparison(results)
    
    print("\n" + "="*100)
    print(f"All visualizations saved to: {analyzer.output_dir}/")
    print("="*100)


if __name__ == '__main__':
    main()
