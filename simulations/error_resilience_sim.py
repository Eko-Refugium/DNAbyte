"""
Cost-Normalized Error Resilience Analysis
All encodings: sequence_length=100bp, normalized cost via mean parameter
Uses full Simulation pipeline: Binarize->Encode->Synthesis->Sequencing->Cluster->Consensus
"""

import os
import sys
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Any
import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.params import Params
from simulations.simulation import Simulation


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


class SimulationBasedResilience:
    """Test encodings with full pipeline and cost normalization"""
    
    def __init__(self, output_dir='./simulations/error_resilience_sim'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = []
    
    def run_analysis(self,
                    encodings: List[str],
                    target_cost_bp: int,
                    error_rates: List[float],
                    base_params: Dict[str, Any],
                    encoding_params: Dict[str, Dict[str, Any]],
                    num_runs: int = 3):
        """Run analysis with cost-normalized encodings using full pipeline"""
        
        print("\n" + "="*100)
        print("COST-NORMALIZED ERROR RESILIENCE - FULL PIPELINE")
        print("="*100)
        print(f"Test data: Bohemian_Rhapsody_Lyrics.txt (1950 bytes = 15,600 bits)")
        print(f"Sequence length: 100 bp (all encodings)")
        print(f"Target cost: {target_cost_bp:,} bp (normalized via mean parameter)")
        print(f"Error rates: {[f'{r*100:.0f}%' for r in error_rates]}")
        print(f"Runs per condition: {num_runs}")
        print(f"Pipeline: Binarize -> Encode -> Synthesis(mean) -> Sequencing(kmere) -> Cluster -> Consensus")
        print("="*100 + "\n")
        
        # Step 1: Calculate baseline costs and mean values
        print("Step 1: Calculating baseline costs and mean values...\n")
        
        mean_values = {}
        baseline_costs = {}
        
        for encoding in encodings:
            try:
                # Run with mean=1 to get baseline cost
                params_dict = base_params.copy()
                params_dict.update(encoding_params[encoding])
                params_dict['encoding_method'] = encoding
                params_dict['sequence_length'] = 100  # All 100bp
                params_dict['mean'] = 1  # Baseline
                params_dict['name'] = f"cost_baseline_{encoding}"
                
                p = Params(**params_dict)
                p.name = params_dict['name']
                
                # Run simulation silently
                sim = Simulation([p])
                results = sim.run()
                
                # Extract cost from results
                baseline_cost = target_cost_bp
                if results:
                    # Results is a dict, get the first (only) simulation result
                    for sim_name, sim_result in results.items():
                        if 'step2' in sim_result and 'number_of_codewords' in sim_result['step2']:
                            # Rough cost estimate: codewords * sequence_length
                            num_codewords = len(p.sequence_length) if hasattr(p, 'sequence_length') else 100
                            baseline_cost = num_codewords * 100  # Approximate
                        break
                
                baseline_costs[encoding] = baseline_cost
                
                # Calculate mean needed to reach target cost
                mean_needed = max(1, int(target_cost_bp / baseline_cost))
                mean_values[encoding] = mean_needed
                
                print(f"  {encoding:18s}: baseline={baseline_cost:6,d}bp  →  mean={mean_needed:3d} for {target_cost_bp:,}bp")
                
            except Exception as e:
                print(f"  {encoding:18s}: ERROR - {str(e)[:60]}")
                mean_values[encoding] = 1
                baseline_costs[encoding] = target_cost_bp
        
        # Step 2: Run error resilience tests
        print("\n" + "="*100)
        print("Step 2: Running error resilience tests with cost-normalized parameters...\n")
        
        for encoding in tqdm(encodings, desc="Encodings"):
            for error_rate in tqdm(error_rates, desc=f"{encoding} ({mean_values[encoding]}x mean)", leave=False):
                successful = 0
                
                for run in range(num_runs):
                    try:
                        # Create simulation with normalized cost
                        params_dict = base_params.copy()
                        params_dict.update(encoding_params[encoding])
                        params_dict['encoding_method'] = encoding
                        params_dict['sequence_length'] = 100  # All 100bp
                        params_dict['mean'] = mean_values[encoding]  # Normalized cost
                        params_dict['kmer_p_sub'] = error_rate  # Substitution error rate
                        params_dict['kmer_p_ins'] = error_rate * 0.5  # Insertion
                        params_dict['kmer_p_del'] = error_rate * 0.5  # Deletion
                        params_dict['name'] = f"{encoding}_err{error_rate:.2f}_run{run}"
                        
                        p = Params(**params_dict)
                        p.name = params_dict['name']
                        
                        # Run simulation
                        sim = Simulation([p])
                        results = sim.run()
                        
                        # Check if recovery was successful
                        if results:
                            for sim_name, sim_result in results.items():
                                if sim_result.get('status') == 'SUCCESS':
                                    successful += 1
                                break
                        
                    except Exception as e:
                        pass  # Count as failure
                
                cost_bp = baseline_costs[encoding] * mean_values[encoding]
                
                result = ResilienceResult(
                    encoding=encoding,
                    error_rate=error_rate,
                    cost_bp=cost_bp,
                    successful_runs=successful,
                    total_runs=num_runs
                )
                
                self.results.append(result)
                
                print(f"  {encoding:18s} @ {error_rate*100:5.1f}% error: "
                      f"{result.recovery_rate*100:5.1f}% ({successful}/{num_runs}), cost={cost_bp:7,d}bp")
        
        # Visualize and summarize
        self._visualize()
        self._print_summary()
    
    def _visualize(self):
        """Create visualization graphs"""
        
        if not self.results:
            print("No results to visualize")
            return
        
        error_rates_set = sorted(set(r.error_rate for r in self.results))
        encodings_set = sorted(set(r.encoding for r in self.results))
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
        
        # Graph 1: Recovery vs Error Rate (all normalized to same cost)
        for enc in encodings_set:
            data = sorted([r for r in self.results if r.encoding == enc],
                         key=lambda x: x.error_rate)
            if data:
                rates = [r.error_rate * 100 for r in data]
                recoveries = [r.recovery_rate * 100 for r in data]
                ax1.plot(rates, recoveries, marker='o', label=enc, linewidth=2.5, markersize=8)
        
        ax1.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax1.set_title('Recovery vs Error Rate (All Equal Cost)', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10, loc='best')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim([-5, 105])
        ax1.set_xlim([-1, 21])
        
        # Graph 2: Recovery by encoding at highest error rate
        if error_rates_set:
            max_error = max(error_rates_set)
            data_at_max = sorted([r for r in self.results if r.error_rate == max_error],
                                key=lambda x: x.recovery_rate, reverse=True)
            
            if data_at_max:
                names = [r.encoding for r in data_at_max]
                recoveries = [r.recovery_rate * 100 for r in data_at_max]
                colors = plt.cm.viridis(np.linspace(0, 1, len(names)))
                
                ax2.bar(names, recoveries, color=colors, alpha=0.8, edgecolor='black', linewidth=1.5)
                ax2.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
                ax2.set_title(f'Recovery at Highest Error Rate ({max_error*100:.0f}%)', 
                             fontsize=14, fontweight='bold')
                ax2.set_ylim([0, 105])
                ax2.grid(True, alpha=0.3, axis='y')
                
                # Add value labels
                for i, (name, rec) in enumerate(zip(names, recoveries)):
                    ax2.text(i, rec + 2, f'{rec:.0f}%', ha='center', fontweight='bold')
        
        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, 'resilience_normalized.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"\n✓ Saved: {plot_file}")
    
    def _print_summary(self):
        """Print results summary"""
        
        print("\n" + "="*100)
        print("RESULTS SUMMARY - COST NORMALIZED")
        print("="*100)
        
        error_rates_set = sorted(set(r.error_rate for r in self.results))
        encodings_set = sorted(set(r.encoding for r in self.results))
        
        for error_rate in error_rates_set:
            print(f"\nError Rate: {error_rate*100:.1f}%")
            data = sorted([r for r in self.results if r.error_rate == error_rate],
                         key=lambda x: x.recovery_rate, reverse=True)
            for r in data:
                print(f"  {r.encoding:18s}: {r.recovery_rate*100:5.1f}% ({r.successful_runs}/{r.total_runs}), cost={r.cost_bp:7,d}bp")


if __name__ == '__main__':
    
    encodings = ['goldman', 'church', 'gcplus',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    # Parameters for 100bp sequences
    encoding_params = {
        'goldman': {
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'church': {
            'rs_num': 10,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'gcplus': {
            'gcplus_k': 8,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'max_density': {
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'no_homopolymer': {
            'max_homopolymer': 4,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'wukong': {
            'rs_num': 10,
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'yinyang': {
            'add_primer': False,
            'primer_length': 0,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
    }
    
    base_params = {
        'filename': 'Bohemian_Rhapsody_Lyrics.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'kmere',
        'kmer_k': 1,
        'kmer_seed': 42,
        'storage_conditions': None,
    }
    
    target_cost_bp = 20000  # All encodings will be normalized to this cost
    error_rates = [0.01, 0.05, 0.10]  # 1%, 5%, 10% error
    
    tester = SimulationBasedResilience()
    tester.run_analysis(
        encodings=encodings,
        target_cost_bp=target_cost_bp,
        error_rates=error_rates,
        base_params=base_params,
        encoding_params=encoding_params,
        num_runs=3  # 3 runs per condition
    )
    
    print("\n✓ Complete!")
