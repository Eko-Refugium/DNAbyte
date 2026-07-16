"""
Simulation: Sequencing error impact on encoding strategies
Demonstrates how different encoding methods respond to varying sequencing error rates.
Shows that encoding strategy performance depends on sequencing technology choice.

This simulation generates data for paper Figure showing sequencing error dependency.
"""
import os
import matplotlib.pyplot as plt
from datetime import datetime
import pickle
import numpy as np

from dnabyte.params import Params
from simulations.simulation import Simulation


def run_sequencing_error_impact():
    """
    Compare encoding strategies across different sequencing error profiles.
    Identical conditions except for sequencing error rates.
    """
    
    encoding_methods = ['max_density', 'church', 'wukong']
    
    # Sequencing error rates to test (realistic range - higher for visible results)
    error_rates = [0.05, 0.10, 0.15, 0.20]  # 5% to 20%
    
    repeats = 3
    
    # Create parameter configurations
    params_list = []
    for error_rate in error_rates:
        for encoding_method in encoding_methods:
            for repeat in range(repeats):
                params = Params(
                    name=f'seq_error_{error_rate}_{encoding_method}_{repeat+1}',
                    filename='textfile_40b.txt',
                    
                    # Vary encoding and sequencing error
                    encoding_method=encoding_method,
                    sequencing_error_rate=error_rate,  # Variable sequencing error
                    
                    # Fixed conditions
                    binarization_method='default',
                    synthesis_method=None,
                    storage_conditions='biogene',
                    sequencing_method='iid',           # IID with mixed error types
                    iid_substitution_rate=error_rate * 0.7,
                    iid_insertion_rate=error_rate * 0.15,
                    iid_deletion_rate=error_rate * 0.15,
                    years=10,
                    
                    codeword_length=500,
                    dna_barcode_length=75,
                )
                params_list.append(params)
    
    print(f"Running {len(params_list)} simulations: {len(encoding_methods)} encodings × {len(error_rates)} error rates × {repeats} repeats")
    
    # Run simulations
    sim = Simulation(params_list)
    results = sim.run(paralel=False)
    
    # Aggregate results by encoding method and error rate
    aggregated = {}
    for encoding_method in encoding_methods:
        aggregated[encoding_method] = {rate: {'SUCCESS': 0, 'FAILURE': 0} for rate in error_rates}
    
    for key, data in results.items():
        # Extract parameters from key
        for error_rate in error_rates:
            for encoding_method in encoding_methods:
                if f'{error_rate}' in key and encoding_method in key:
                    status = data.get('status', 'UNKNOWN')
                    if status in aggregated[encoding_method][error_rate]:
                        aggregated[encoding_method][error_rate][status] += 1
    
    # Calculate success rates
    success_data = {}
    for encoding_method in encoding_methods:
        success_data[encoding_method] = []
        for error_rate in error_rates:
            total = aggregated[encoding_method][error_rate]['SUCCESS'] + aggregated[encoding_method][error_rate]['FAILURE']
            if total > 0:
                rate = aggregated[encoding_method][error_rate]['SUCCESS'] / total
            else:
                rate = 0
            success_data[encoding_method].append(rate)
    
    print("\n" + "="*70)
    print("SEQUENCING ERROR IMPACT ON ENCODING STRATEGIES")
    print("="*70)
    for encoding_method in encoding_methods:
        print(f"\n{encoding_method}:")
        for error_rate, success_rate in zip(error_rates, success_data[encoding_method]):
            print(f"  Error rate {error_rate:.3f}: {success_rate:.1%} success")
    
    # Generate visualization
    fig, ax = plt.subplots(figsize=(11, 7))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    markers = ['o', 's', '^']
    
    for (encoding_method, color, marker) in zip(encoding_methods, colors, markers):
        ax.plot([e*100 for e in error_rates], success_data[encoding_method], 
                marker=marker, linewidth=2.5, markersize=8, label=encoding_method, 
                color=color, alpha=0.8)
    
    ax.set_xlabel('Sequencing Error Rate (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reconstruction Success Rate', fontsize=12, fontweight='bold')
    ax.set_title('Sequencing Error Impact on Encoding Strategy Performance\n(Identical Synthesis, Storage, and Recovery)', 
                 fontsize=14, fontweight='bold')
    ax.set_ylim([-0.05, 1.05])
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='best', framealpha=0.9)
    
    # Save figure
    job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_path = os.path.join('simulations', 'simlogs', f'sequencing_error_impact_{job_identifier}.png')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")
    
    # Save results
    pickle_path = os.path.join('simulations', 'simlogs', f'sequencing_error_impact_{job_identifier}.pickle')
    with open(pickle_path, 'wb') as f:
        pickle.dump({'results': results, 'aggregated': aggregated, 'success_data': success_data}, f)
    print(f"Results saved: {pickle_path}")
    
    return results, aggregated, success_data


if __name__ == '__main__':
    results, aggregated, success_data = run_sequencing_error_impact()
