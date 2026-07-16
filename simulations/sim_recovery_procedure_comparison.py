"""
Simulation: Recovery procedure comparison
Demonstrates how different clustering and recovery methods affect reconstruction performance.
Shows that recovery method choice influences end-to-end system performance.

This simulation generates data for paper Figure showing recovery procedure comparison.
"""
import os
import matplotlib.pyplot as plt
from datetime import datetime
import pickle
import numpy as np

from dnabyte.params import Params
from simulations.simulation import Simulation


def run_recovery_procedure_comparison():
    """
    Compare different recovery procedures (clustering + recovery combinations).
    Tests how recovery method choice affects reconstruction success.
    """
    
    encoding_method = 'church'  # Fixed encoding
    
    # Recovery procedures (clustering + recovery combinations)
    recovery_procedures = [
        ('primer_grouper', 'debruijn'),    # Primer-based clustering + De Bruijn recovery
        # Add more recovery procedures here as they're implemented
    ]
    
    # Sequencing error rates to test recovery robustness
    error_rates = [0, 0.05, 0.10, 0.15]  # 5%, 10%, 15% error rates
    
    repeats = 3
    
    # Create parameter configurations
    params_list = []
    for (clustering_method, recovery_method) in recovery_procedures:
        for error_rate in error_rates:
            for repeat in range(repeats):
                params = Params(
                    name=f'recovery_{clustering_method}_{recovery_method}_{error_rate}_{repeat+1}',
                    filename='textfile_40b.txt',
                    
                    # Fixed encoding with synthesis to create copies for recovery
                    encoding_method=encoding_method,
                    synthesis_method='mesa',           # MESA synthesis creates multiple copies
                    
                    # Vary recovery procedure and sequencing error
                    clustering_method=clustering_method,
                    recovery_method=recovery_method,
                    sequencing_error_rate=error_rate,
                    
                    # Fixed conditions
                    binarization_method='default',
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
    
    print(f"Running {len(params_list)} simulations: {len(recovery_procedures)} recovery procedures × {len(error_rates)} error rates × {repeats} repeats")
    
    # Run simulations
    sim = Simulation(params_list)
    results = sim.run(paralel=False)
    
    # Aggregate results
    aggregated = {}
    for (clustering_method, recovery_method) in recovery_procedures:
        procedure_name = f'{clustering_method}+{recovery_method}'
        aggregated[procedure_name] = {rate: {'SUCCESS': 0, 'FAILURE': 0} for rate in error_rates}
    
    for key, data in results.items():
        for (clustering_method, recovery_method) in recovery_procedures:
            for error_rate in error_rates:
                if clustering_method in key and recovery_method in key and f'{error_rate}' in key:
                    procedure_name = f'{clustering_method}+{recovery_method}'
                    status = data.get('status', 'UNKNOWN')
                    if status in aggregated[procedure_name][error_rate]:
                        aggregated[procedure_name][error_rate][status] += 1
    
    # Calculate success rates
    success_data = {}
    for (clustering_method, recovery_method) in recovery_procedures:
        procedure_name = f'{clustering_method}+{recovery_method}'
        success_data[procedure_name] = []
        for error_rate in error_rates:
            total = aggregated[procedure_name][error_rate]['SUCCESS'] + aggregated[procedure_name][error_rate]['FAILURE']
            if total > 0:
                rate = aggregated[procedure_name][error_rate]['SUCCESS'] / total
            else:
                rate = 0
            success_data[procedure_name].append(rate)
    
    print("\n" + "="*70)
    print("RECOVERY PROCEDURE COMPARISON")
    print("="*70)
    for procedure_name in success_data.keys():
        print(f"\n{procedure_name}:")
        for error_rate, success_rate in zip(error_rates, success_data[procedure_name]):
            print(f"  Error rate {error_rate:.3f}: {success_rate:.1%} success")
    
    # Generate visualization
    fig, ax = plt.subplots(figsize=(11, 7))
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    markers = ['o', 's', '^', 'D']
    
    for (procedure_name, color, marker) in zip(success_data.keys(), colors, markers):
        ax.plot([e*100 for e in error_rates], success_data[procedure_name],
                marker=marker, linewidth=2.5, markersize=8, label=procedure_name,
                color=color, alpha=0.8)
    
    ax.set_xlabel('Sequencing Error Rate (%)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Reconstruction Success Rate', fontsize=12, fontweight='bold')
    ax.set_title(f'Recovery Procedure Performance Comparison\n({encoding_method} encoding, identical synthesis and storage)',
                 fontsize=14, fontweight='bold')
    ax.set_ylim([-0.05, 1.05])
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(fontsize=11, loc='best', framealpha=0.9)
    
    # Save figure
    job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_path = os.path.join('simulations', 'simlogs', f'recovery_comparison_{job_identifier}.png')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved: {output_path}")
    
    # Save results
    pickle_path = os.path.join('simulations', 'simlogs', f'recovery_comparison_{job_identifier}.pickle')
    with open(pickle_path, 'wb') as f:
        pickle.dump({'results': results, 'aggregated': aggregated, 'success_data': success_data}, f)
    print(f"Results saved: {pickle_path}")
    
    return results, aggregated, success_data


if __name__ == '__main__':
    results, aggregated, success_data = run_recovery_procedure_comparison()
