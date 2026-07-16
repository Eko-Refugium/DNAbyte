"""
Simulation: Comparative evaluation of encoding strategies
Demonstrates how different encoding methods perform under identical conditions,
showing that encoding strategy performance depends on pipeline configuration.

This simulation generates data for paper Figure showing encoding strategy comparison.
"""
import os
import matplotlib.pyplot as plt
from datetime import datetime
import pickle
import numpy as np

from dnabyte.params import Params
from simulations.simulation import Simulation


def run_encoding_comparison():
    """
    Compare multiple encoding strategies under identical synthesis/storage/sequencing conditions.
    Only the encoding method varies while all other pipeline components remain constant.
    """
    
    # Define encoding methods to compare
    encoding_methods = ['max_density', 'church', 'wukong']
    
    # Number of repeats for statistical robustness
    repeats = 3
    
    # Sequencing error rates for realistic stress testing (higher for visible differences)
    error_rates = [0.05, 0.10, 0.15]  # 5%, 10%, 15% error rates
    
    # Create parameter configurations: one per encoding method, repeated
    params_list = []
    for encoding_method in encoding_methods:
        for error_rate in error_rates:
            for repeat in range(repeats):
                params = Params(
                    name=f'encoding_comparison_{encoding_method}_err{error_rate}_{repeat+1}',
                    filename='textfile_40b.txt',
                    
                    # Vary only encoding method - all else constant
                    encoding_method=encoding_method,
                    binarization_method='default',
                    sequencing_error_rate=error_rate,  # REALISTIC error conditions
                    
                    # Identical synthesis, storage, sequencing conditions
                    synthesis_method=None,
                    storage_conditions='biogene',      # Standard storage
                    sequencing_method='iid',           # IID with substitutions + indels
                    iid_substitution_rate=error_rate * 0.7,  # 70% substitutions
                    iid_insertion_rate=error_rate * 0.15,    # 15% insertions
                    iid_deletion_rate=error_rate * 0.15,     # 15% deletions
                    years=10,                          # Short-term storage
                    
                    # Fixed parameters for reproducibility
                    codeword_length=500,
                    dna_barcode_length=75,
                )
                params_list.append(params)
    
    print(f"Running {len(params_list)} simulations: {len(encoding_methods)} encodings × {len(error_rates)} error rates × {repeats} repeats")
    
    # Run simulations
    sim = Simulation(params_list)
    results = sim.run(paralel=False)
    
    # Aggregate results by encoding method
    encoding_results = {method: {rate: {'SUCCESS': 0, 'FAILURE': 0} for rate in error_rates} for method in encoding_methods}
    
    for key, data in results.items():
        # Extract encoding method and error rate from key
        for method in encoding_methods:
            for error_rate in error_rates:
                if method in key and f'err{error_rate}' in key:
                    status = data.get('status', 'UNKNOWN')
                    if status in encoding_results[method][error_rate]:
                        encoding_results[method][error_rate][status] += 1
                    break
    
    # Calculate success rates
    success_rates = {method: [] for method in encoding_methods}
    for method in encoding_methods:
        for error_rate in sorted(error_rates):
            total = encoding_results[method][error_rate]['SUCCESS'] + encoding_results[method][error_rate]['FAILURE']
            if total > 0:
                rate = encoding_results[method][error_rate]['SUCCESS'] / total
            else:
                rate = 0
            success_rates[method].append(rate)
    
    print("\n" + "="*70)
    print("ENCODING STRATEGY COMPARISON RESULTS")
    print("="*70)
    for method in encoding_methods:
        print(f"\n{method}:")
        for error_rate, success_rate in zip(sorted(error_rates), success_rates[method]):
            print(f"  Error rate {error_rate:.1%}: {success_rate:.1%} success")
    
    # Generate visualization
    fig, ax = plt.subplots(figsize=(11, 7))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    markers = ['o', 's', '^']
    
    for (method, color, marker) in zip(encoding_methods, colors, markers):
        ax.plot([e*100 for e in sorted(error_rates)], success_rates[method],
                marker=marker, linewidth=2.5, markersize=8, label=method,
                color=color, alpha=0.8)
    
    ax.set_xlabel('Sequencing Error Rate (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reconstruction Success Rate', fontsize=12, fontweight='bold')
    ax.set_title('Encoding Strategy Performance Under Sequencing Errors\n(Identical Synthesis, Storage, and Recovery)', 
                 fontsize=14, fontweight='bold')
    ax.set_ylim([-0.05, 1.05])
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='best', framealpha=0.9)
    
    # Save figure
    job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_path = os.path.join('simulations', 'simlogs', f'encoding_comparison_{job_identifier}.png')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")
    
    # Save results
    pickle_path = os.path.join('simulations', 'simlogs', f'encoding_comparison_{job_identifier}.pickle')
    with open(pickle_path, 'wb') as f:
        pickle.dump({'results': results, 'aggregated': encoding_results, 'success_rates': success_rates}, f)
    print(f"Results saved: {pickle_path}")
    
    return results, encoding_results, success_rates


if __name__ == '__main__':
    results, encoding_results, success_rates = run_encoding_comparison()
