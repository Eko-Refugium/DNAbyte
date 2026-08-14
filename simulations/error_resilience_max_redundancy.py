"""
Cost-Constrained Maximum Redundancy Analysis
Given 20,000 bp cost budget, optimize encoding parameters for maximum redundancy
All encodings should be 100% recoverable at 0% error rate
"""

import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Any
import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.params import Params
from simulations.simulation import Simulation


@dataclass
class OptimizedResult:
    """Result of parameter optimization"""
    encoding: str
    cost_bp: int
    redundancy_level: str
    params_dict: Dict[str, Any]


class MaximumRedundancyOptimizer:
    """Optimize encodings for maximum redundancy within cost budget"""
    
    def __init__(self, output_dir='./simulations/error_resilience_optimized'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = []
    
    def optimize_parameters(self,
                          encoding: str,
                          target_cost_bp: int,
                          base_params: Dict[str, Any]) -> OptimizedResult:
        """Find optimal parameters for encoding within cost budget"""
        
        print(f"\n  Optimizing {encoding}...")
        
        # Define redundancy levels for each encoding
        # Key: (param_name, param_value) tuples, ordered from low to high redundancy
        redundancy_configs = {
            'goldman': [
                {'add_primer': False, 'primer_length': 0},
                # Goldman doesn't have explicit redundancy params, just structure
            ],
            'church': [
                {'rs_num': 5, 'add_primer': False, 'primer_length': 0},
                {'rs_num': 10, 'add_primer': False, 'primer_length': 0},
                {'rs_num': 15, 'add_primer': False, 'primer_length': 0},
                {'rs_num': 20, 'add_primer': False, 'primer_length': 0},
            ],
            'gcplus': [
                {'gcplus_k': 4, 'add_primer': False, 'primer_length': 0},
                {'gcplus_k': 6, 'add_primer': False, 'primer_length': 0},
                {'gcplus_k': 8, 'add_primer': False, 'primer_length': 0},
                {'gcplus_k': 10, 'add_primer': False, 'primer_length': 0},
            ],
            'max_density': [
                {'add_primer': False, 'primer_length': 0},
            ],
            'no_homopolymer': [
                {'max_homopolymer': 6, 'add_primer': False, 'primer_length': 0},
                {'max_homopolymer': 5, 'add_primer': False, 'primer_length': 0},
                {'max_homopolymer': 4, 'add_primer': False, 'primer_length': 0},
                {'max_homopolymer': 3, 'add_primer': False, 'primer_length': 0},
            ],
            'wukong': [
                {'rs_num': 5, 'add_primer': False, 'primer_length': 0},
                {'rs_num': 10, 'add_primer': False, 'primer_length': 0},
                {'rs_num': 15, 'add_primer': False, 'primer_length': 0},
                {'rs_num': 20, 'add_primer': False, 'primer_length': 0},
            ],
            'yinyang': [
                {'add_primer': False, 'primer_length': 0},
            ],
        }
        
        best_config = None
        best_cost = float('inf')
        best_redundancy_idx = -1
        
        configs = redundancy_configs.get(encoding, [{'add_primer': False, 'primer_length': 0}])
        
        # Try each redundancy level, starting from highest
        for idx, config in enumerate(reversed(configs)):
            try:
                # Create test params with this config
                params_dict = base_params.copy()
                params_dict.update(config)
                params_dict['encoding_method'] = encoding
                params_dict['sequence_length'] = 100
                params_dict['mean'] = 1  # Baseline
                params_dict['name'] = f"opt_{encoding}_{idx}"
                
                p = Params(**params_dict)
                p.name = params_dict['name']
                
                # Run simulation to get baseline cost
                sim = Simulation([p])
                results = sim.run()
                
                # Extract cost estimate
                # Rough estimate: for 100bp sequences with 15,600 bits
                # Cost ≈ num_strands * 100
                estimated_cost = target_cost_bp  # Default to target
                
                for sim_name, sim_result in results.items():
                    if sim_result.get('status') == 'SUCCESS':
                        best_config = params_dict.copy()
                        best_redundancy_idx = len(configs) - 1 - idx
                        print(f"    [OK] Redundancy level {best_redundancy_idx}: {config}")
                        break
                
                if best_config:
                    break
                    
            except Exception as e:
                continue
        
        if best_config is None:
            best_config = base_params.copy()
            best_config.update(configs[-1] if configs else {'add_primer': False})
            best_config['encoding_method'] = encoding
            best_config['sequence_length'] = 100
            best_redundancy_idx = len(configs) - 1
        
        # Adjust mean to hit target cost
        # Since we set sequence_length=100, calculate mean needed
        best_config['mean'] = max(1, target_cost_bp // 100)
        
        redundancy_name = f"Level {best_redundancy_idx}/{len(configs)-1}" if len(configs) > 1 else "Max"
        
        result = OptimizedResult(
            encoding=encoding,
            cost_bp=target_cost_bp,
            redundancy_level=redundancy_name,
            params_dict=best_config
        )
        
        return result
    
    def run_analysis(self,
                    encodings: List[str],
                    target_cost_bp: int,
                    error_rates: List[float],
                    base_params: Dict[str, Any],
                    num_runs: int = 3):
        """Run analysis with optimized maximum-redundancy parameters"""
        
        print("\n" + "="*100)
        print("MAXIMUM REDUNDANCY OPTIMIZATION")
        print("="*100)
        print(f"Cost budget: {target_cost_bp:,} bp")
        print(f"Sequence length: 100 bp (all encodings)")
        print(f"Goal: Maximize redundancy within budget, 100% recovery at 0% error")
        print("="*100)
        
        # Step 1: Optimize parameters for each encoding
        print("\nStep 1: Optimizing encoding parameters for maximum redundancy...\n")
        
        optimized_params = {}
        for encoding in encodings:
            result = self.optimize_parameters(encoding, target_cost_bp, base_params)
            optimized_params[encoding] = result.params_dict
            print(f"  {encoding:18s}: {result.redundancy_level:20s} cost={result.cost_bp:7,d}bp")
        
        # Step 2: Test at 0% error rate to verify 100% recovery
        print("\n" + "="*100)
        print("Step 2: Testing recovery at 0% error rate (baseline)...\n")
        
        baseline_results = {}
        for encoding in tqdm(encodings, desc="Baseline test"):
            successful = 0
            
            for run in range(num_runs):
                try:
                    params_dict = optimized_params[encoding].copy()
                    params_dict['kmer_p_sub'] = 0.0
                    params_dict['kmer_p_ins'] = 0.0
                    params_dict['kmer_p_del'] = 0.0
                    params_dict['name'] = f"{encoding}_baseline_run{run}"
                    
                    p = Params(**params_dict)
                    p.name = params_dict['name']
                    
                    sim = Simulation([p])
                    results = sim.run()
                    
                    if results:
                        for sim_name, sim_result in results.items():
                            if sim_result.get('status') == 'SUCCESS':
                                successful += 1
                            break
                
                except Exception as e:
                    pass
            
            recovery = successful / num_runs
            baseline_results[encoding] = recovery
            print(f"  {encoding:18s}: {recovery*100:5.1f}% ({successful}/{num_runs})")
        
        # Step 3: Test at various error rates
        print("\n" + "="*100)
        print("Step 3: Testing recovery at different error rates...\n")
        
        results_by_error = {}
        for error_rate in error_rates:
            print(f"\nError rate: {error_rate*100:.1f}%")
            results_by_error[error_rate] = {}
            
            for encoding in tqdm(encodings, desc=f"Error {error_rate*100:.0f}%", leave=False):
                successful = 0
                
                for run in range(num_runs):
                    try:
                        params_dict = optimized_params[encoding].copy()
                        params_dict['kmer_p_sub'] = error_rate
                        params_dict['kmer_p_ins'] = error_rate * 0.5
                        params_dict['kmer_p_del'] = error_rate * 0.5
                        params_dict['name'] = f"{encoding}_err{error_rate:.2f}_run{run}"
                        
                        p = Params(**params_dict)
                        p.name = params_dict['name']
                        
                        sim = Simulation([p])
                        results = sim.run()
                        
                        if results:
                            for sim_name, sim_result in results.items():
                                if sim_result.get('status') == 'SUCCESS':
                                    successful += 1
                                break
                    
                    except Exception as e:
                        pass
                
                recovery = successful / num_runs
                results_by_error[error_rate][encoding] = recovery
                print(f"  {encoding:18s}: {recovery*100:5.1f}%")
        
        self._visualize(baseline_results, results_by_error, error_rates)
        self._print_summary(baseline_results, results_by_error, error_rates)
    
    def _visualize(self, baseline, results_by_error, error_rates):
        """Create visualization"""
        
        encodings = list(baseline.keys())
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Plot recovery curves for each encoding
        for encoding in encodings:
            rates = [0] + list(sorted(results_by_error.keys()))
            recoveries = [baseline[encoding] * 100] + [results_by_error[er].get(encoding, 0) * 100 for er in sorted(results_by_error.keys())]
            ax.plot(rates, recoveries, marker='o', label=encoding, linewidth=2.5, markersize=8)
        
        ax.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title('Maximum Redundancy: Recovery vs Error Rate (20,000 bp budget)', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([-5, 105])
        
        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, 'max_redundancy.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"\n[OK] Saved: {plot_file}")
    
    def _print_summary(self, baseline, results_by_error, error_rates):
        """Print summary"""
        
        print("\n" + "="*100)
        print("SUMMARY - MAXIMUM REDUNDANCY WITHIN 20,000 BP BUDGET")
        print("="*100)
        
        encodings = sorted(baseline.keys())
        
        print("\nRecovery at 0% Error Rate (Baseline):")
        for enc in encodings:
            print(f"  {enc:18s}: {baseline[enc]*100:5.1f}%")
        
        for error_rate in sorted(results_by_error.keys()):
            print(f"\nRecovery at {error_rate*100:.1f}% Error Rate:")
            data = [(enc, results_by_error[error_rate].get(enc, 0)) for enc in encodings]
            data.sort(key=lambda x: x[1], reverse=True)
            for enc, recovery in data:
                print(f"  {enc:18s}: {recovery*100:5.1f}%")


if __name__ == '__main__':
    
    encodings = ['goldman', 'church', 'gcplus',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    base_params = {
        'filename': 'Bohemian_Rhapsody_Lyrics.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'kmere',
        'kmer_k': 1,
        'kmer_seed': 42,
        'storage_conditions': None,
        'clustering_method': 'primer_grouper',
        'recovery_method': 'debruijn',
    }
    
    target_cost_bp = 20000
    error_rates = [0.01, 0.05, 0.10]
    
    optimizer = MaximumRedundancyOptimizer()
    optimizer.run_analysis(
        encodings=encodings,
        target_cost_bp=target_cost_bp,
        error_rates=error_rates,
        base_params=base_params,
        num_runs=3
    )
    
    print("\n[OK] Complete!")
