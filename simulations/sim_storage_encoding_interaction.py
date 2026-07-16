"""
Simulation: Storage duration and encoding strategy interaction
Demonstrates that encoding strategy performance depends on storage conditions and duration.
Shows how long-term storage affects different encoding methods differently.

This simulation generates data for paper Figure showing storage-encoding interaction.
"""
import os
import matplotlib.pyplot as plt
from datetime import datetime
import pickle
import numpy as np

from dnabyte.params import Params
from simulations.simulation import Simulation


def run_storage_encoding_interaction():
    """
    Compare encoding strategies across different storage durations.
    Shows how storage-induced degradation affects each encoding differently.
    """
    
    encoding_methods = ['max_density', 'church', 'wukong']
    
    # Storage durations (years) - logarithmic scale
    storage_durations = [1, 10, 100, 1000, 10000]
    
    repeats = 4
    
    # Create parameter configurations
    params_list = []
    for storage_duration in storage_durations:
        for encoding_method in encoding_methods:
            for repeat in range(repeats):
                params = Params(
                    name=f'storage_encoding_{storage_duration}yr_{encoding_method}_{repeat+1}',
                    filename='textfile_40b.txt',
                    
                    # Vary encoding and storage duration
                    encoding_method=encoding_method,
                    years=storage_duration,  # Variable storage duration
                    sequencing_error_rate=0.05,  # REALISTIC sequencing errors (5%)
                    
                    # Fixed conditions
                    binarization_method='default',
                    synthesis_method=None,
                    storage_conditions='biogene',      # Realistic degradation model
                    sequencing_method='iid',           # IID with mixed error types
                    iid_substitution_rate=0.035,       # 3.5% substitutions
                    iid_insertion_rate=0.0075,         # 0.75% insertions
                    iid_deletion_rate=0.0075,          # 0.75% deletions
                    
                    codeword_length=500,
                    dna_barcode_length=75,
                )
                params_list.append(params)
    
    print(f"Running {len(params_list)} simulations: {len(encoding_methods)} encodings × {len(storage_durations)} durations × {repeats} repeats")
    
    # Run simulations
    sim = Simulation(params_list)
    results = sim.run(paralel=False)
    
    # Aggregate results
    aggregated = {}
    for encoding_method in encoding_methods:
        aggregated[encoding_method] = {duration: {'SUCCESS': 0, 'FAILURE': 0} for duration in storage_durations}
    
    for key, data in results.items():
        for duration in storage_durations:
            for encoding_method in encoding_methods:
                if f'{duration}yr' in key and encoding_method in key:
                    status = data.get('status', 'UNKNOWN')
                    if status in aggregated[encoding_method][duration]:
                        aggregated[encoding_method][duration][status] += 1
    
    # Calculate success rates
    success_data = {}
    for encoding_method in encoding_methods:
        success_data[encoding_method] = []
        for duration in storage_durations:
            total = aggregated[encoding_method][duration]['SUCCESS'] + aggregated[encoding_method][duration]['FAILURE']
            if total > 0:
                rate = aggregated[encoding_method][duration]['SUCCESS'] / total
            else:
                rate = 0
            success_data[encoding_method].append(rate)
    
    print("\n" + "="*70)
    print("STORAGE DURATION vs ENCODING STRATEGY")
    print("="*70)
    for encoding_method in encoding_methods:
        print(f"\n{encoding_method}:")
        for duration, success_rate in zip(storage_durations, success_data[encoding_method]):
            print(f"  {duration:5d} years: {success_rate:.1%} success")
    
    # Generate visualization
    fig, ax = plt.subplots(figsize=(11, 7))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    markers = ['o', 's', '^']
    
    for (encoding_method, color, marker) in zip(encoding_methods, colors, markers):
        ax.plot(storage_durations, success_data[encoding_method],
                marker=marker, linewidth=2.5, markersize=8, label=encoding_method,
                color=color, alpha=0.8)
    
    ax.set_xlabel('Storage Duration (years)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reconstruction Success Rate', fontsize=12, fontweight='bold')
    ax.set_title('Storage Duration Impact on Encoding Strategy Performance\n(Biogene Storage Conditions, Perfect Sequencing)',
                 fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.set_ylim([-0.05, 1.05])
    ax.grid(True, alpha=0.3, linestyle='--', which='both')
    ax.legend(fontsize=11, loc='best', framealpha=0.9)
    
    # Save figure
    job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_path = os.path.join('simulations', 'simlogs', f'storage_encoding_interaction_{job_identifier}.png')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")
    
    # Save results
    pickle_path = os.path.join('simulations', 'simlogs', f'storage_encoding_interaction_{job_identifier}.pickle')
    with open(pickle_path, 'wb') as f:
        pickle.dump({'results': results, 'aggregated': aggregated, 'success_data': success_data}, f)
    print(f"Results saved: {pickle_path}")
    
    return results, aggregated, success_data


if __name__ == '__main__':
    results, aggregated, success_data = run_storage_encoding_interaction()
