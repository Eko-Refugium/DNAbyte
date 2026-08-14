"""
Cost-Normalized Error Resilience Analysis using Full Pipeline
All encodings at sequence_length=200, normalized to same cost via mean parameter
Full pipeline: Encode > Synthesize (mean copies) > Sequencing (kmere errors) > Cluster+Consensus > Decode
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
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.sequence import SimulateSequencing
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus


@dataclass
class ResilienceResult:
    """Cost-normalized resilience result"""
    encoding: str
    error_rate: float
    target_cost_bp: int
    mean_copies: int
    successful_runs: int
    total_runs: int
    actual_cost_bp: int
    
    @property
    def recovery_rate(self) -> float:
        return self.successful_runs / self.total_runs
    
    @property
    def effective_density(self) -> float:
        """bits per nucleotide at given error rate and cost"""
        return (15600 / max(self.actual_cost_bp, 1)) * self.recovery_rate


class CostNormalizedResilience:
    """Test encodings with cost equalization via synthesis mean"""
    
    def __init__(self, output_dir='./simulations/error_analysis_normalized'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = []
    
    def _calculate_base_cost(self, encoding: str, params_dict: Dict[str, Any]) -> int:
        """Calculate base DNA cost at mean=1"""
        try:
            params = Params(**params_dict)
            test_data = Data(file_paths=['./tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'])
            binary_code = Binarize(params).binarize(test_data)
            dna_codewords, _ = Encode(params).encode(binary_code)
            
            cost_bp = sum(len(seq) for seq in dna_codewords.data)
            return cost_bp
        except Exception as e:
            return 0
    
    def _calculate_mean_for_target_cost(self, base_cost: int, target_cost: int) -> int:
        """Calculate mean (synthesis copies) needed to reach target cost"""
        if base_cost == 0:
            return 1
        mean = max(1, int(np.ceil(target_cost / base_cost)))
        return mean
    
    def _run_full_pipeline(self,
                          encoding: str,
                          error_rate: float,
                          mean: int,
                          params_dict: Dict[str, Any]) -> bool:
        """
        Run full pipeline: Encode > Synthesize > Sequencing > Cluster+Consensus
        Returns: True if recovery succeeded, False otherwise
        """
        try:
            # Add error and synthesis params
            params = Params(
                encoding_method=encoding,
                binarization_method='default',
                sequence_length=200,
                synthesis_method='nosynthpoly',
                mean=mean,
                std_dev=0,
                sequencing_method='kmere' if error_rate > 0 else None,
                kmer_k=1,
                kmer_p_ins=error_rate * 0.5,
                kmer_p_del=error_rate * 0.25,
                kmer_p_sub=error_rate * 0.25,
                kmer_seed=42,
                clustering_method='primer_grouper',
                recovery_method='debruijn',
                **{k: v for k, v in params_dict.items() 
                   if k not in ['encoding_method', 'binarization_method', 'sequence_length',
                                'synthesis_method', 'mean', 'std_dev', 'sequencing_method',
                                'clustering_method', 'recovery_method']}
            )
            
            # Step 1: Binarize
            test_data = Data(file_paths=['./tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'])
            binary_code = Binarize(params).binarize(test_data)
            
            # Step 2: Encode
            enc = Encode(params)
            dna_codewords, _ = enc.encode(binary_code)
            
            # Step 3: Synthesize (create copies)
            dna_obj = InSilicoDNA(dna_codewords.data)
            syn = SimulateSynthesis(params)
            dna_syn, _ = syn.simulate(dna_obj)
            
            # Step 4: Sequencing with errors
            if error_rate > 0:
                seq = SimulateSequencing(params)
                dna_seq, _ = seq.simulate(dna_syn)
            else:
                dna_seq = dna_syn
            
            # Step 5: Recovery via Cluster+Consensus
            cluster_obj = Cluster(params)
            clustered_data, _ = cluster_obj.cluster(dna_seq)
            
            consensus_obj = Consensus(params)
            recovered_data, _ = consensus_obj.call(clustered_data)
            
            # Step 6: Decode
            data_dec, valid, _ = enc.decode(recovered_data)
            
            return valid
            
        except Exception as e:
            return False
    
    def test_encoding(self,
                     encoding: str,
                     error_rate: float,
                     target_cost_bp: int,
                     base_params: Dict[str, Any],
                     working_params: Dict[str, Any],
                     num_runs: int = 20) -> ResilienceResult:
        """Test encoding with cost normalization"""
        
        # Build params
        params_dict = base_params.copy()
        if encoding in working_params:
            params_dict.update(working_params[encoding])
        params_dict['encoding_method'] = encoding
        params_dict['filename'] = 'Bohemian_Rhapsody_Lyrics.txt'
        
        # Calculate base cost
        base_cost = self._calculate_base_cost(encoding, params_dict)
        if base_cost == 0:
            return None
        
        # Calculate mean for target cost
        mean = self._calculate_mean_for_target_cost(base_cost, target_cost_bp)
        actual_cost = base_cost * mean
        
        # Run tests
        successful = 0
        for run in range(num_runs):
            if self._run_full_pipeline(encoding, error_rate, mean, params_dict):
                successful += 1
        
        result = ResilienceResult(
            encoding=encoding,
            error_rate=error_rate,
            target_cost_bp=target_cost_bp,
            mean_copies=mean,
            successful_runs=successful,
            total_runs=num_runs,
            actual_cost_bp=actual_cost
        )
        
        self.results.append(result)
        return result
    
    def run_analysis(self,
                    encodings: List[str],
                    target_cost_bp: int,
                    error_rates: List[float],
                    base_params: Dict[str, Any],
                    working_params: Dict[str, Any],
                    num_runs: int = 20):
        """Run full cost-normalized analysis"""
        
        print("\n" + "="*120)
        print("COST-NORMALIZED ERROR RESILIENCE ANALYSIS")
        print("="*120)
        print(f"Target Cost: {target_cost_bp:,} bp (same for all encodings via mean parameter)")
        print(f"Error Rates: {[f'{r*100:.1f}%' for r in error_rates]}")
        print(f"Runs per condition: {num_runs}")
        print(f"Pipeline: Encode(seq_len=200) > Synthesize(mean copies) > Sequencing(kmere) > Cluster+Consensus > Decode")
        print("="*120 + "\n")
        
        # First pass: calculate base costs and mean factors
        print("STEP 1: CALCULATING BASE COSTS AND MEAN FACTORS\n")
        base_costs = {}
        mean_factors = {}
        
        for encoding in encodings:
            params_dict = base_params.copy()
            if encoding in working_params:
                params_dict.update(working_params[encoding])
            params_dict['encoding_method'] = encoding
            params_dict['filename'] = 'Bohemian_Rhapsody_Lyrics.txt'
            
            base_cost = self._calculate_base_cost(encoding, params_dict)
            mean = self._calculate_mean_for_target_cost(base_cost, target_cost_bp)
            actual_cost = base_cost * mean
            base_costs[encoding] = base_cost
            mean_factors[encoding] = mean
            
            print(f"  {encoding:18s}: base={base_cost:6d}bp, mean={mean:2d}, actual={actual_cost:7d}bp")
        
        print("\nSTEP 2: TESTING ERROR RESILIENCE\n")
        
        # Run tests
        for encoding in tqdm(encodings, desc="Encodings"):
            for error_rate in tqdm(error_rates, desc=encoding, leave=False):
                result = self.test_encoding(
                    encoding, error_rate, target_cost_bp,
                    base_params, working_params, num_runs
                )
                
                if result:
                    print(f"{result.encoding:18s} @ {result.error_rate*100:5.1f}%: "
                          f"mean={result.mean_copies:2d}, {result.recovery_rate*100:5.1f}% recovery "
                          f"({result.successful_runs:2d}/{result.total_runs}), "
                          f"cost={result.actual_cost_bp:7d}bp")
        
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
        
        # Graph 1: Recovery vs Error Rate (cost-normalized)
        for enc in encodings_set:
            data = sorted([r for r in self.results if r.encoding == enc],
                         key=lambda x: x.error_rate)
            if data:
                rates = [r.error_rate * 100 for r in data]
                recoveries = [r.recovery_rate * 100 for r in data]
                ax1.plot(rates, recoveries, marker='o', label=enc, linewidth=2.5, markersize=8)
        
        ax1.set_xlabel('Sequencing Error Rate (%)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Recovery Rate (%)', fontsize=12, fontweight='bold')
        ax1.set_title('Cost-Normalized Recovery (All at ~30kb via mean)', fontsize=13, fontweight='bold')
        ax1.legend(fontsize=10, loc='best')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim([-5, 105])
        
        # Graph 2: Effective Density vs Error Rate
        for enc in encodings_set:
            data = sorted([r for r in self.results if r.encoding == enc],
                         key=lambda x: x.error_rate)
            if data:
                rates = [r.error_rate * 100 for r in data]
                density = [r.effective_density for r in data]
                ax2.plot(rates, density, marker='s', label=enc, linewidth=2.5, markersize=8)
        
        ax2.set_xlabel('Sequencing Error Rate (%)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Effective Density (bits/nt)', fontsize=12, fontweight='bold')
        ax2.set_title('Data Density vs Error Rate (Cost-Normalized)', fontsize=13, fontweight='bold')
        ax2.legend(fontsize=10, loc='best')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plot_file = os.path.join(self.output_dir, 'cost_normalized_analysis.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"\n[OK] Saved: {plot_file}")
    
    def _print_summary(self):
        """Print results summary"""
        
        print("\n" + "="*120)
        print("COST-NORMALIZED RESULTS SUMMARY")
        print("="*120)
        
        error_rates_set = sorted(set(r.error_rate for r in self.results))
        encodings_set = sorted(set(r.encoding for r in self.results))
        
        for error_rate in error_rates_set:
            print(f"\n{error_rate*100:.1f}% Sequencing Error Rate:")
            data = sorted([r for r in self.results if r.error_rate == error_rate],
                         key=lambda x: x.recovery_rate, reverse=True)
            for r in data:
                print(f"  {r.encoding:18s}: mean={r.mean_copies:2d}, {r.recovery_rate*100:5.1f}% recovery "
                      f"({r.successful_runs:2d}/{r.total_runs}), cost={r.actual_cost_bp:7d}bp, "
                      f"eff_dens={r.effective_density:.4f}")
        
        print("\n" + "="*120)
        print("BEST ENCODING BY ERROR RATE")
        print("="*120)
        
        for error_rate in error_rates_set:
            data = sorted([r for r in self.results if r.error_rate == error_rate],
                         key=lambda x: x.recovery_rate, reverse=True)
            if data:
                best = data[0]
                print(f"{error_rate*100:5.1f}%: {best.encoding:18s} ({best.recovery_rate*100:5.1f}%)")


if __name__ == '__main__':
    
    encodings = ['goldman', 'church', 'gcplus',
                 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    # All fixed at sequence_length=200
    working_params = {
        'goldman': {
            'sequence_length': 200,
            'add_primer': False,
            'primer_length': 0,
        },
        'church': {
            'sequence_length': 200,
            'rs_num': 10,
            'add_primer': False,
            'primer_length': 0,
        },
        'gcplus': {
            'sequence_length': 200,
            'gcplus_k': 10,
            'add_primer': False,
            'primer_length': 0,
        },
        'max_density': {
            'codeword_length': 200,
            'add_primer': False,
            'primer_length': 0,
        },
        'no_homopolymer': {
            'codeword_length': 200,
            'max_homopolymer': 5,
            'add_primer': False,
            'primer_length': 0,
        },
        'wukong': {
            'sequence_length': 200,
            'rs_num': 10,
            'add_primer': False,
            'primer_length': 0,
        },
        'yinyang': {
            'sequence_length': 200,
            'add_primer': False,
            'primer_length': 0,
        },
    }
    
    base_params = {
        'encoding_method': None,
        'binarization_method': 'default',
    }
    
    # Target cost: 30,000 bp (all encodings normalized to this)
    target_cost_bp = 30000
    
    error_rates = [0.01, 0.05, 0.10, 0.15, 0.20]
    
    analyzer = CostNormalizedResilience()
    analyzer.run_analysis(
        encodings=encodings,
        target_cost_bp=target_cost_bp,
        error_rates=error_rates,
        base_params=base_params,
        working_params=working_params,
        num_runs=20
    )
    
    print("\n[OK] Complete!")
