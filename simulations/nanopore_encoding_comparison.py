"""
Nanopore Sequencing Error Resilience Comparison
Compare different encodings' resilience to nanopore sequencing errors
Track insertions, deletions, and substitutions for each encoding
"""

import os
import sys
import json
import random
import string
from dataclasses import dataclass, asdict
from typing import Dict, List, Any, Tuple
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from dnabyte.params import Params
from simulations.simulation import Simulation


@dataclass
class ErrorStats:
    """Error statistics for a single run"""
    encoding: str
    insertions: int
    deletions: int
    substitutions: int
    total_bases: int
    recovery_success: bool


def create_random_file(output_path: str, size_bytes: int = 500):
    """Create a random text file for encoding"""
    content = ''.join(random.choices(string.ascii_letters + string.digits + ' \n', k=size_bytes))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(content)
    return output_path


def align_sequences(seq1: str, seq2: str) -> Tuple[int, int, int]:
    """
    Simple alignment to count insertions, deletions, and substitutions
    Uses Needleman-Wunsch algorithm for global alignment
    Returns: (insertions, deletions, substitutions)
    """
    # Simple DP alignment
    m, n = len(seq1), len(seq2)
    
    # Create DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j],    # deletion
                    dp[i][j-1],    # insertion
                    dp[i-1][j-1]   # substitution
                )
    
    # Backtrack to count error types
    insertions = 0
    deletions = 0
    substitutions = 0
    
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and seq1[i-1] == seq2[j-1]:
            i -= 1
            j -= 1
        elif i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + 1:
            # substitution
            substitutions += 1
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + 1:
            # deletion
            deletions += 1
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j-1] + 1:
            # insertion
            insertions += 1
            j -= 1
        else:
            break
    
    return insertions, deletions, substitutions


def count_sequence_errors(original_sequences: List[str], modified_sequences: List[str]) -> Tuple[int, int, int, int]:
    """
    Count total insertions, deletions, and substitutions across all sequences
    Returns: (insertions, deletions, substitutions, total_bases)
    """
    total_insertions = 0
    total_deletions = 0
    total_substitutions = 0
    total_bases = 0
    
    # Flatten if nested
    def flatten(lst):
        result = []
        for item in lst:
            if isinstance(item, list):
                result.extend(flatten(item))
            else:
                result.append(item)
        return result
    
    orig_flat = flatten(original_sequences)
    mod_flat = flatten(modified_sequences)
    
    # Align each pair of sequences
    for orig, mod in zip(orig_flat, mod_flat):
        if isinstance(orig, str) and isinstance(mod, str):
            ins, dels, subs = align_sequences(orig, mod)
            total_insertions += ins
            total_deletions += dels
            total_substitutions += subs
            total_bases += len(orig)
    
    return total_insertions, total_deletions, total_substitutions, total_bases


class NanoporeEncodingComparison:
    """Compare encoding resilience to nanopore sequencing errors"""
    
    def __init__(self, output_dir='./simulations/nanopore_comparison'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.results = []
    
    def run_comparison(self,
                      encodings: List[str],
                      nanopore_methods: List[int] = [39, 40],  # 39=1D, 40=2D
                      num_runs: int = 10):
        """
        Run comparison for all encodings with both nanopore methods
        """
        
        print("\n" + "="*100)
        print("NANOPORE ENCODING COMPARISON")
        print("="*100)
        
        # Create random test file
        test_file_path = './tests/testfiles/random_test_nanopore.txt'
        create_random_file(test_file_path, size_bytes=500)
        print(f"\nCreated random test file: {test_file_path}")
        
        # Define encoding parameters
        encoding_params = self._get_encoding_params()
        
        # Base parameters
        base_params = {
            'filename': 'random_test_nanopore.txt',
            'binarization_method': 'default',
            'synthesis_method': 'nosynthpoly',
            'sequencing_method': 'mesa',
            'storage_conditions': None,
            'sequence_length': 100,
            'mean': 10,
        }
        
        # Store all results for plotting
        all_results = []
        
        # Test each encoding with each nanopore method
        for encoding in encodings:
            print(f"\n{'='*80}")
            print(f"Testing encoding: {encoding.upper()}")
            print(f"{'='*80}")
            
            for nanopore_id in nanopore_methods:
                method_name = "Nanopore 1D" if nanopore_id == 39 else "Nanopore 2D"
                print(f"\n  Sequencing method: {method_name}")
                
                # Run multiple times for statistics
                run_results = []
                for run_idx in range(num_runs):
                    try:
                        # Create params for this run
                        test_params = base_params.copy()
                        test_params.update(encoding_params[encoding])
                        test_params['encoding_method'] = encoding
                        test_params['mesa_sequencing_id'] = nanopore_id
                        test_params['name'] = f"{encoding}_nanopore_{nanopore_id}_run_{run_idx}"
                        test_params['seed'] = 42 + run_idx
                        
                        # Run simulation with custom error tracking
                        insertions, deletions, substitutions, total_bases, recovery_success = self._run_single_simulation(test_params)
                        
                        run_results.append({
                            'encoding': encoding,
                            'nanopore_method': method_name,
                            'run': run_idx,
                            'insertions': insertions,
                            'deletions': deletions,
                            'substitutions': substitutions,
                            'total_bases': total_bases,
                            'recovery_success': recovery_success,
                            'status': 'SUCCESS'
                        })
                        
                        print(f"    Run {run_idx+1}/{num_runs}: Ins={insertions}, Del={deletions}, Sub={substitutions}, Success={recovery_success}")
                                
                    except Exception as e:
                        print(f"    Run {run_idx} failed: {str(e)}")
                        run_results.append({
                            'encoding': encoding,
                            'nanopore_method': method_name,
                            'run': run_idx,
                            'insertions': 0,
                            'deletions': 0,
                            'substitutions': 0,
                            'total_bases': 0,
                            'recovery_success': False,
                            'status': 'FAILURE'
                        })
                
                # Calculate success rate and average errors for this encoding + method combo
                success_count = sum(1 for r in run_results if r.get('recovery_success', False))
                success_rate = success_count / num_runs if num_runs > 0 else 0
                
                # Calculate average error rates
                avg_ins = np.mean([r.get('insertions', 0) for r in run_results])
                avg_del = np.mean([r.get('deletions', 0) for r in run_results])
                avg_sub = np.mean([r.get('substitutions', 0) for r in run_results])
                avg_bases = np.mean([r.get('total_bases', 1) for r in run_results])
                
                print(f"    Success rate: {success_rate*100:.1f}% ({success_count}/{num_runs})")
                print(f"    Average errors: Ins={avg_ins:.1f}, Del={avg_del:.1f}, Sub={avg_sub:.1f} (out of {avg_bases:.0f} bases)")
                
                all_results.extend(run_results)
        
        # Save results
        results_file = os.path.join(self.output_dir, 'comparison_results.json')
        with open(results_file, 'w') as f:
            json.dump(all_results, f, indent=2)
        
        print(f"\n\nResults saved to: {results_file}")
        
        # Print summary table
        self.print_summary_table(all_results)
        
        # Plot results
        self.plot_results(all_results)
        
        return all_results
    
    def print_summary_table(self, results: List[Dict]):
        """Print a summary table of results"""
        print("\n" + "="*100)
        print("SUMMARY TABLE")
        print("="*100)
        
        # Group by encoding
        encodings = sorted(set(r['encoding'] for r in results))
        
        print(f"\n{'Encoding':<15} {'Method':<15} {'Success Rate':<15} {'Insertions %':<13} {'Deletions %':<13} {'Substitutions %':<15}")
        print("-" * 100)
        
        for encoding in encodings:
            for method_name in ['Nanopore 1D', 'Nanopore 2D']:
                enc_results = [r for r in results if r['encoding'] == encoding and r['nanopore_method'] == method_name]
                
                if enc_results:
                    success_rate = np.mean([r.get('recovery_success', False) for r in enc_results]) * 100
                    # Calculate error percentages (errors / total_bases * 100)
                    avg_ins_pct = np.mean([r.get('insertions', 0) / max(r.get('total_bases', 1), 1) * 100 for r in enc_results])
                    avg_del_pct = np.mean([r.get('deletions', 0) / max(r.get('total_bases', 1), 1) * 100 for r in enc_results])
                    avg_sub_pct = np.mean([r.get('substitutions', 0) / max(r.get('total_bases', 1), 1) * 100 for r in enc_results])
                    
                    print(f"{encoding:<15} {method_name:<15} {success_rate:>6.1f}%        {avg_ins_pct:>8.2f}%    {avg_del_pct:>8.2f}%    {avg_sub_pct:>8.2f}%")
        
        print("="*100)
    
    def _run_single_simulation(self, test_params: Dict) -> Tuple[int, int, int, int, bool]:
        """
        Run a single simulation and extract error statistics from Mesa simulator
        Returns: (insertions, deletions, substitutions, total_bases, recovery_success)
        """
        from dnabyte.data_classes.base import Data
        from dnabyte.data_classes.insilicodna import InSilicoDNA
        from dnabyte.data_classes.nucleobasecode import NucleobaseCode
        from dnabyte.binarize import Binarize
        from dnabyte.encode import Encode
        from dnabyte.synthesize import SimulateSynthesis
        from dnabyte.sequence import SimulateSequencing
        
        # Create params
        p = Params(**test_params)
        p.name = test_params['name']
        
        # Run encoding pipeline manually to capture sequences
        try:
            # Step 1: Binarize
            bin_obj = Binarize(p)
            file_paths = ['./tests/testfiles/' + p.filename]
            data_obj = Data(file_paths=file_paths)
            binary_code = bin_obj.binarize(data_obj)
            
            # Step 2: Encode
            encoder = Encode(p)
            nucleobase_code, _ = encoder.encode(binary_code)
            
            # Step 3: Synthesis
            # Ensure nucleobase_code is the right type
            if not isinstance(nucleobase_code, NucleobaseCode):
                # Convert InSilicoDNA to NucleobaseCode if needed
                if isinstance(nucleobase_code, InSilicoDNA):
                    nucleobase_code = NucleobaseCode(nucleobase_code.data)
            
            synthesizer = SimulateSynthesis(p)
            synthesized_dna, _ = synthesizer.simulate(nucleobase_code)
            
            # Calculate total bases before sequencing
            original_sequences = self._extract_sequences(synthesized_dna)
            total_bases = sum(len(seq) for seq in original_sequences if isinstance(seq, str))
            
            # Step 4: Sequencing (where errors are introduced)
            # Mesa simulator now returns error counts in the info dict
            sequencer = SimulateSequencing(p)
            sequenced_dna, seq_info = sequencer.simulate(synthesized_dna)
            
            # Extract error counts from Mesa simulator info
            insertions = seq_info.get('total_insertions', 0)
            deletions = seq_info.get('total_deletions', 0)
            substitutions = seq_info.get('total_substitutions', 0)
            
            # Continue with decoding to check recovery
            # Step 5: Process
            processed_code, _ = encoder.process(sequenced_dna)
            
            # Step 6: Decode
            decoded_data, valid, _ = encoder.decode(processed_code)
            
            # Step 7: Compare
            comparison, _ = decoded_data.compare(decoded_data, binary_code)
            recovery_success = (comparison != 'ERROR')
            
            return insertions, deletions, substitutions, total_bases, recovery_success
            
        except Exception as e:
            print(f"      Error in simulation: {str(e)}")
            import traceback
            traceback.print_exc()
            return 0, 0, 0, 0, False
    
    def _extract_sequences(self, dna_object) -> List[str]:
        """Extract DNA sequences from InSilicoDNA object"""
        sequences = []
        
        if hasattr(dna_object, 'data'):
            data = dna_object.data
            
            # Handle different data structures
            if isinstance(data, list):
                for item in data:
                    if isinstance(item, str):
                        sequences.append(item)
                    elif isinstance(item, list):
                        sequences.extend(self._extract_sequences_from_list(item))
                    elif hasattr(item, 'sequence'):
                        sequences.append(item.sequence)
            elif isinstance(data, str):
                sequences.append(data)
        
        return sequences
    
    def _extract_sequences_from_list(self, data_list) -> List[str]:
        """Recursively extract sequences from nested list"""
        sequences = []
        for item in data_list:
            if isinstance(item, str):
                sequences.append(item)
            elif isinstance(item, list):
                sequences.extend(self._extract_sequences_from_list(item))
            elif hasattr(item, 'sequence'):
                sequences.append(item.sequence)
        return sequences
    
    def _get_encoding_params(self) -> Dict[str, Dict[str, Any]]:
        """Get parameters for each encoding"""
        return {
            'goldman': {
                'add_primer': True,
                'primer_length': 20,
                'clustering_method': 'kmere_cluster',
                'recovery_method': 'debruijn',
                'mean': 5
            },
            'church': {
                'rs_num': 10,
                'max_homopolymer': 4,
                'add_redundancy': True,
                'add_primer': True,
                'primer_length': 20,
                'clustering_method': 'kmere_cluster',
                'recovery_method': 'debruijn',
                'mean': 5
            },
            'gcplus': {
                'gcplus_k': 100,
                'gcplus_l': 10,
                'gcplus_c1': 4,
                'add_primer': True,
                'primer_length': 20,
                'clustering_method': 'kmere_cluster',
                'recovery_method': 'debruijn',
                'mean': 5
            },
            'max_density': {
                'dna_barcode_length': 20,
                'codeword_maxlength_positions': 20,
                'codeword_length': 200,
                'mean': 5,
                'clustering_method': 'kmere_cluster',
                'recovery_method': 'debruijn'
            },
            'no_homopolymer': {
                'dna_barcode_length': 20,
                'codeword_maxlength_positions': 20,
                'codeword_length': 200,
                'mean': 5,
                'clustering_method': 'kmere_cluster',
                'recovery_method': 'debruijn'
            },
            'wukong': {
                'rs_num': 3,
                'max_homopolymer': 3,
                'min_gc': 0.45,
                'max_gc': 0.55,
                'rule_num': 1,
                'add_redundancy': True,
                'add_primer': True,
                'primer_length': 20,
                'clustering_method': 'kmere_cluster',
                'recovery_method': 'debruijn',
                'mean': 5
            },
            'yinyang': {
                'max_homopolymer': 4,
                'yinyang_search_count': 1000,
                'mean': 5,                
                'clustering_method': 'kmere_cluster',
                'recovery_method': 'debruijn'
            },
        }
    
    def plot_results(self, results: List[Dict]):
        """Plot comparison bar charts for both success rates and error types"""
        
        # Aggregate results by encoding and nanopore method
        success_data = defaultdict(lambda: {'1D': [], '2D': []})
        error_data = defaultdict(lambda: {
            '1D': {'insertions': [], 'deletions': [], 'substitutions': []},
            '2D': {'insertions': [], 'deletions': [], 'substitutions': []}
        })
        
        for result in results:
            encoding = result['encoding']
            method = result['nanopore_method']
            success = 1 if result.get('recovery_success', False) else 0
            
            # Get error rates (as percentages)
            total_bases = result.get('total_bases', 1)
            if total_bases == 0:
                total_bases = 1
            
            ins_rate = (result.get('insertions', 0) / total_bases) * 100
            del_rate = (result.get('deletions', 0) / total_bases) * 100
            sub_rate = (result.get('substitutions', 0) / total_bases) * 100
            
            if 'Nanopore 1D' in method:
                success_data[encoding]['1D'].append(success)
                error_data[encoding]['1D']['insertions'].append(ins_rate)
                error_data[encoding]['1D']['deletions'].append(del_rate)
                error_data[encoding]['1D']['substitutions'].append(sub_rate)
            elif 'Nanopore 2D' in method:
                success_data[encoding]['2D'].append(success)
                error_data[encoding]['2D']['insertions'].append(ins_rate)
                error_data[encoding]['2D']['deletions'].append(del_rate)
                error_data[encoding]['2D']['substitutions'].append(sub_rate)
        
        encodings = sorted(success_data.keys())
        
        # Plot 1: Success Rates
        fig, ax = plt.subplots(figsize=(12, 6))
        
        success_1d = [np.mean(success_data[enc]['1D']) * 100 if success_data[enc]['1D'] else 0 for enc in encodings]
        success_2d = [np.mean(success_data[enc]['2D']) * 100 if success_data[enc]['2D'] else 0 for enc in encodings]
        
        x = np.arange(len(encodings))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, success_1d, width, label='Nanopore 1D', alpha=0.8, color='#2ecc71')
        bars2 = ax.bar(x + width/2, success_2d, width, label='Nanopore 2D', alpha=0.8, color='#3498db')
        
        ax.set_xlabel('Encoding Method', fontsize=12, fontweight='bold')
        ax.set_ylabel('Recovery Success Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title('Encoding Resilience to Nanopore Sequencing Errors', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(encodings, rotation=45, ha='right')
        ax.legend()
        ax.set_ylim(0, 105)
        ax.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.1f}%',
                       ha='center', va='bottom', fontsize=8)
        
        plt.tight_layout()
        plot_path = os.path.join(self.output_dir, 'nanopore_success_rates.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"\nSuccess rate plot saved to: {plot_path}")
        plt.close()
        
        # Plot 2: Error Types (Combined 1D and 2D average)
        fig, ax = plt.subplots(figsize=(14, 7))
        
        # Calculate average error rates across both methods
        avg_insertions = []
        avg_deletions = []
        avg_substitutions = []
        
        for enc in encodings:
            all_ins = error_data[enc]['1D']['insertions'] + error_data[enc]['2D']['insertions']
            all_del = error_data[enc]['1D']['deletions'] + error_data[enc]['2D']['deletions']
            all_sub = error_data[enc]['1D']['substitutions'] + error_data[enc]['2D']['substitutions']
            
            avg_insertions.append(np.mean(all_ins) if all_ins else 0)
            avg_deletions.append(np.mean(all_del) if all_del else 0)
            avg_substitutions.append(np.mean(all_sub) if all_sub else 0)
        
        x = np.arange(len(encodings))
        width = 0.25
        
        bars1 = ax.bar(x - width, avg_insertions, width, label='Insertions', alpha=0.8, color='#e74c3c')
        bars2 = ax.bar(x, avg_deletions, width, label='Deletions', alpha=0.8, color='#f39c12')
        bars3 = ax.bar(x + width, avg_substitutions, width, label='Substitutions', alpha=0.8, color='#9b59b6')
        
        ax.set_xlabel('Encoding Method', fontsize=12, fontweight='bold')
        ax.set_ylabel('Error Rate (%)', fontsize=12, fontweight='bold')
        ax.set_title('Nanopore Sequencing Error Types by Encoding', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(encodings, rotation=45, ha='right')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bars in [bars1, bars2, bars3]:
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2., height,
                           f'{height:.2f}%',
                           ha='center', va='bottom', fontsize=7)
        
        plt.tight_layout()
        plot_path = os.path.join(self.output_dir, 'nanopore_error_types.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"Error types plot saved to: {plot_path}")
        plt.close()
        
        # Plot 3: Detailed comparison for 1D vs 2D
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))
        
        # Nanopore 1D errors
        ins_1d = [np.mean(error_data[enc]['1D']['insertions']) if error_data[enc]['1D']['insertions'] else 0 for enc in encodings]
        del_1d = [np.mean(error_data[enc]['1D']['deletions']) if error_data[enc]['1D']['deletions'] else 0 for enc in encodings]
        sub_1d = [np.mean(error_data[enc]['1D']['substitutions']) if error_data[enc]['1D']['substitutions'] else 0 for enc in encodings]
        
        x = np.arange(len(encodings))
        width = 0.25
        
        ax1.bar(x - width, ins_1d, width, label='Insertions', alpha=0.8, color='#e74c3c')
        ax1.bar(x, del_1d, width, label='Deletions', alpha=0.8, color='#f39c12')
        ax1.bar(x + width, sub_1d, width, label='Substitutions', alpha=0.8, color='#9b59b6')
        
        ax1.set_xlabel('Encoding Method', fontsize=11, fontweight='bold')
        ax1.set_ylabel('Error Rate (%)', fontsize=11, fontweight='bold')
        ax1.set_title('Nanopore 1D Sequencing Errors', fontsize=12, fontweight='bold')
        ax1.set_xticks(x)
        ax1.set_xticklabels(encodings, rotation=45, ha='right')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)
        
        # Nanopore 2D errors
        ins_2d = [np.mean(error_data[enc]['2D']['insertions']) if error_data[enc]['2D']['insertions'] else 0 for enc in encodings]
        del_2d = [np.mean(error_data[enc]['2D']['deletions']) if error_data[enc]['2D']['deletions'] else 0 for enc in encodings]
        sub_2d = [np.mean(error_data[enc]['2D']['substitutions']) if error_data[enc]['2D']['substitutions'] else 0 for enc in encodings]
        
        ax2.bar(x - width, ins_2d, width, label='Insertions', alpha=0.8, color='#e74c3c')
        ax2.bar(x, del_2d, width, label='Deletions', alpha=0.8, color='#f39c12')
        ax2.bar(x + width, sub_2d, width, label='Substitutions', alpha=0.8, color='#9b59b6')
        
        ax2.set_xlabel('Encoding Method', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Error Rate (%)', fontsize=11, fontweight='bold')
        ax2.set_title('Nanopore 2D Sequencing Errors', fontsize=12, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(encodings, rotation=45, ha='right')
        ax2.legend()
        ax2.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plot_path = os.path.join(self.output_dir, 'nanopore_detailed_comparison.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"Detailed comparison plot saved to: {plot_path}")
        plt.close()


if __name__ == '__main__':
    
    # Encodings to test (you can comment out some to test faster)
    encodings = [
        'goldman',
        'church',
        'gcplus',
        'max_density',
        'no_homopolymer',
        'wukong',
        'yinyang'
    ]
    
    # Nanopore methods: 39 = 1D, 40 = 2D
    nanopore_methods = [39, 40]
    
    # Number of runs per encoding+method combination
    # Start with 3-5 for quick testing, increase to 10-20 for final results
    num_runs = 5
    
    print("\n" + "="*100)
    print("NANOPORE ENCODING COMPARISON SIMULATION")
    print("="*100)
    print(f"\nEncodings to test: {len(encodings)}")
    print(f"Nanopore methods: {len(nanopore_methods)} (1D and 2D)")
    print(f"Runs per combination: {num_runs}")
    print(f"Total simulations: {len(encodings) * len(nanopore_methods) * num_runs}")
    print("\nThis may take several minutes to complete...")
    print("="*100)
    
    # Run comparison
    comparator = NanoporeEncodingComparison()
    results = comparator.run_comparison(
        encodings=encodings,
        nanopore_methods=nanopore_methods,
        num_runs=num_runs
    )
    
    print("\n" + "="*100)
    print("[COMPLETE] Nanopore encoding comparison finished!")
    print(f"Results and plots saved in: {comparator.output_dir}")
    print("="*100)
