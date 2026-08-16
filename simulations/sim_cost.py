"""
Cost-Constrained Maximum Redundancy Analysis (Simplified)
All encodings set to 20,000 bp with robust pre-configured parameters
Test recovery at 0% and higher error rates
"""

import os
import sys
from dataclasses import dataclass
from typing import Dict, List, Any
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.params import Params
from simulations.simulation import Simulation


@dataclass
class RecoveryResult:
    """Recovery test result"""
    encoding: str
    error_rate: float
    recovery_rate: float
    successful_runs: int
    total_runs: int


class RobustParameterAnalyzer:
    """Test encodings with robust pre-configured parameters"""
    
    def __init__(self, output_dir='./simulations/error_resilience_robust'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = []

    def vary_parameters(self, base_params: Dict[str, Any], encoding_params: Dict[str, Dict[str, Any]]):
        """Search for parameters with 100% recovery at 0% error, scaled to 20,000 bp target"""
        
        optimized_params = {}
        target_cost_scaled = 1000  # At mean=10, target is 20x smaller: 20,000 * (10/200) = 1,000 bp
        
        for encoding, base_config in encoding_params.items():
            print(f"\n  Optimizing {encoding}...")
            
            working_configs = []  # Configs that achieve 100% recovery
            
            # Incrementally increase redundancy level
            for redundancy_level in range(0, 100, 5):
                config = self._create_config_at_level(encoding, redundancy_level)
                
                test_params = base_params.copy()
                test_params.update(config)
                test_params['encoding_method'] = encoding
                test_params['sequence_length'] = 100
                test_params['mean'] = 10  # Small mean for fast testing
                test_params['kmer_p_sub'] = 0.0
                test_params['kmer_p_ins'] = 0.0
                test_params['kmer_p_del'] = 0.0
                test_params['name'] = f"{encoding}_opt_test"
                
                
                try:
                    p = Params(**test_params)
                    p.name = test_params['name']
                    sim = Simulation([p])
                    results = sim.run()
                    
                    bitstreams_match = False
                    cost_bp = 0
                    if results:
                        for sim_name, sim_result in results.items():
                            # Check if bitstreams match from step9 comparison
                            step9 = sim_result.get('step9', {})
                            if isinstance(step9, dict):
                                bitstreams_match = step9.get('bitstreams_match', False)
                            cost_bp = sim_result.get('total_dna_length', 0)
                            break
                    
                    if bitstreams_match and cost_bp > 0:
                        cost_pct = (cost_bp / target_cost_scaled) * 100
                        print(f"    Level {redundancy_level:2d}: Cost {cost_bp:,} bp ({cost_pct:.1f}%) - WORKS")
                        working_configs.append((config, cost_bp, redundancy_level))
                    else:
                        print(f"    Level {redundancy_level:2d}: FAILED")
                
                except Exception as e:
                    print(f"    Level {redundancy_level:2d}: ERROR - {str(e)[:100]}")
            
            # Among working configs, pick closest to target
            if working_configs:
                working_configs.sort(key=lambda x: abs(x[1] - target_cost_scaled))
                best_config, best_cost, best_level = working_configs[0]
                optimized_params[encoding] = best_config
                cost_pct = (best_cost / target_cost_scaled) * 100
                print(f"    [SELECTED] Level {best_level}: {best_cost:,} bp ({cost_pct:.1f}%)")
            else:
                # Fallback to base config
                optimized_params[encoding] = base_config.copy()
                print(f"    [FALLBACK] Using base config (no working configs found)")
        print(optimized_params)
        breakpoint()
        return optimized_params
    
    def _create_config_at_level(self, encoding: str, level: int) -> Dict[str, Any]:
        """Create parameter config for redundancy level 0-100"""
        
        if encoding == 'goldman':
            return {'add_primer': True, 'primer_length': 20}
        
        elif encoding == 'church':
            rs_num = int(5 + (level / 100) * 15)
            max_homo = max(4, int(6 - (level / 100) * 2))
            add_red = level > 50
            return {
                'rs_num': rs_num,
                'max_homopolymer': max_homo,
                'add_redundancy': add_red,
                'add_primer': False,
                'primer_length': 20
            }
        
        elif encoding == 'gcplus':
            gcplus_k = max(50, int(200 - (level / 100) * 150))
            gcplus_l = min(12, int(8 + (level / 100) * 4))
            gcplus_c1 = 2 + (level // 50)
            return {
                'gcplus_k': gcplus_k,
                'gcplus_l': gcplus_l,
                'gcplus_c1': gcplus_c1,
                'add_primer': False,
                'primer_length': 20
            }
        
        elif encoding == 'max_density':
            maxlen = min(28, int(10 + (level / 100) * 18))
            return {
                'dna_barcode_length': 20,
                'codeword_maxlength_positions': maxlen,
                'codeword_length': 100
            }
        
        elif encoding == 'no_homopolymer':
            maxlen = min(28, int(10 + (level / 100) * 18))
            return {
                'dna_barcode_length': 20,
                'codeword_maxlength_positions': maxlen,
                'codeword_length': 100
            }
        
        elif encoding == 'wukong':
            rs_num = int(5 + (level / 100) * 15)
            max_homo = max(2, int(5 - (level / 100) * 3))
            gc_width = max(0.10, 0.20 - (level / 100) * 0.10)
            gc_center = 0.50
            add_red = level > 40
            return {
                'rs_num': rs_num,
                'max_homopolymer': max_homo,
                'min_gc': gc_center - gc_width/2,
                'max_gc': gc_center + gc_width/2,
                'rule_num': 1 + (level // 50),
                'add_redundancy': add_red,
                'add_primer': False,
                'primer_length': 20
            }
        
        elif encoding == 'yinyang':
            max_homo = max(2, int(6 - (level / 100) * 4))
            search_count = int(50 + (level / 100) * 150)
            return {
                'max_homopolymer': max_homo,
                'yinyang_search_count': search_count,
                'add_primer': False,
                'primer_length': 20
            }
        
        return {'add_primer': True, 'primer_length': 20}
    
    def run_analysis(self,
                    encodings: List[str],
                    target_cost_bp: int,
                    error_rates: List[float],
                    base_params: Dict[str, Any],
                    encoding_params: Dict[str, Dict[str, Any]],
                    num_runs: int = 3):
        """Run analysis with robust parameters optimized for recovery"""
        
        print("\n" + "="*100)
        print("ROBUST PARAMETER ANALYSIS")
        print("="*100)
        print(f"Cost budget: {target_cost_bp:,} bp (all encodings equal)")
        print(f"Sequence length: 100 bp (all encodings)")
        print(f"Parameters: Optimized via search")
        print("="*100)
        
        # Step 0: Optimize parameters for each encoding
        print("\nStep 0: Searching for optimal parameters within budget...\n")
        optimized_params = self.vary_parameters(base_params, encoding_params)
        encoding_params = optimized_params
        
        # Step 1: Test at 0% error (baseline)
        print("\nStep 1: Baseline test at 0% error rate (no errors)...\n")
        
        baseline_results = {}
        for encoding in tqdm(encodings, desc="Baseline"):
            successful = 0
            
            for run in range(num_runs):
                try:
                    params_dict = base_params.copy()
                    params_dict.update(encoding_params[encoding])
                    params_dict['encoding_method'] = encoding
                    params_dict['sequence_length'] = 100
                    params_dict['mean'] = 10  # Normalize to cost
                    params_dict['kmer_p_sub'] = 0.0
                    params_dict['kmer_p_ins'] = 0.0
                    params_dict['kmer_p_del'] = 0.0
                    params_dict['name'] = f"{encoding}_baseline_{run}"
                    
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
            result = RecoveryResult(encoding, 0.0, recovery, successful, num_runs)
            self.results.append(result)
            print(f"  {encoding:18s}: {recovery*100:5.1f}% ({successful}/{num_runs})")
        
        # Step 2: Test at error rates
        print("\n" + "="*100)
        print("Step 2: Testing at various error rates...\n")
        
        for error_rate in sorted(error_rates):
            print(f"Error rate: {error_rate*100:.1f}%")
            
            for encoding in tqdm(encodings, desc=f"Error {error_rate*100:.0f}%", leave=False):
                successful = 0
                
                for run in range(num_runs):
                    try:
                        params_dict = base_params.copy()
                        params_dict.update(encoding_params[encoding])
                        params_dict['encoding_method'] = encoding
                        params_dict['sequence_length'] = 100
                        params_dict['mean'] = 10
                        params_dict['kmer_p_sub'] = error_rate
                        params_dict['kmer_p_ins'] = error_rate * 0.5
                        params_dict['kmer_p_del'] = error_rate * 0.5
                        params_dict['name'] = f"{encoding}_err{error_rate:.2f}_{run}"
                        
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
                result = RecoveryResult(encoding, error_rate, recovery, successful, num_runs)
                self.results.append(result)
                print(f"  {encoding:18s}: {recovery*100:5.1f}%")
        
        self._visualize()
        self._print_summary()
    
    def _visualize(self):
        """Create visualization"""
        
        encodings = sorted(set(r.encoding for r in self.results))
        error_rates = sorted(set(r.error_rate for r in self.results))
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        for encoding in encodings:
            data = sorted([r for r in self.results if r.encoding == encoding],
                         key=lambda x: x.error_rate)
            if data:
                rates = [r.error_rate * 100 for r in data]
                recoveries = [r.recovery_rate * 100 for r in data]
                ax.plot(rates, recoveries, marker='o', label=encoding, linewidth=2.5, markersize=8)
        
        ax.set_xlabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title('Robust Parameters: Recovery vs Error Rate (20,000 bp budget)', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([-5, 105])
        
        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, 'robust_recovery.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"\n[OK] Saved: {plot_file}")
    
    def _print_summary(self):
        """Print summary"""
        
        print("\n" + "="*100)
        print("SUMMARY - RECOVERY WITH ROBUST PARAMETERS")
        print("="*100)
        
        error_rates = sorted(set(r.error_rate for r in self.results))
        
        for error_rate in error_rates:
            print(f"\n{error_rate*100:.1f}% Error Rate:")
            data = sorted([r for r in self.results if r.error_rate == error_rate],
                         key=lambda x: x.recovery_rate, reverse=True)
            for r in data:
                print(f"  {r.encoding:18s}: {r.recovery_rate*100:5.1f}% ({r.successful_runs}/{r.total_runs})")


if __name__ == '__main__':
    
    encodings = ['goldman', 'church', 'gcplus',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    # Robust parameters for each encoding
    encoding_params = {
        'goldman': {
            'add_primer': True,
            'primer_length': 20,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'church': {
            'rs_num': 20,  # High redundancy
            'add_primer': True,
            'primer_length': 20,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'gcplus': {
            'gcplus_k': 4,  # Lower k = higher redundancy
            'add_primer': True,
            'primer_length': 20,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'max_density': {
            'dna_barcode_length': 20,
            'add_primer': True,
            'primer_length': 20,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'no_homopolymer': {
            'max_homopolymer': 3,  # Tighter constraint = more redundancy
            'dna_barcode_length': 20,
            'add_primer': True,
            'primer_length': 20,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'wukong': {
            'rs_num': 20,  # High redundancy
            'add_primer': True,
            'primer_length': 20,
            'clustering_method': 'primer_grouper',
            'recovery_method': 'debruijn'
        },
        'yinyang': {
            'add_primer': True,
            'primer_length': 20,
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
        'recovery_method': 'debruijn',
        'clustering_method': 'kmere_cluster',
        'storage_conditions': None,
        'synthesis_method': 'nosynthpoly',
        'kmer_seed': 42,
        'storage_conditions': None,
    }
    
    target_cost_bp = 20000
    error_rates = [0.0, 0.01, 0.05, 0.10]
    
    analyzer = RobustParameterAnalyzer()
    analyzer.run_analysis(
        encodings=encodings,
        target_cost_bp=target_cost_bp,
        error_rates=error_rates,
        base_params=base_params,
        encoding_params=encoding_params,
        num_runs=3
    )
    
    print("\n[OK] Complete!")


