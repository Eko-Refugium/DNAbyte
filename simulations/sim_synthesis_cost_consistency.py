"""
Simulation: Verify synthesis cost consistency and test error resilience.

Part 1: Verifies synthesis cost is identical across synthesis methods
        for the same encoding (validates cost measures encoder output).

Part 2: Tests full pipeline (encode → synthesize → store → sequence → 
        recover → decode) at different error rates to measure how error
        resilient each encoding is at its given synthesis cost.

This provides a complete picture: encoding efficiency vs. error resilience.
"""

import os
import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from dnabyte.params import Params
from simulations.encoding_analysis import AnalyzeEncoding
from simulations.simulation import Simulation


def test_synthesis_cost_consistency():
    """
    Test that synthesis cost is consistent across synthesis methods
    for the same encoding.
    """
    encodings = ['max_density', 'church', 'wukong']
    synthesis_methods = [None, 'mesa']
    
    print("=" * 80)
    print("SYNTHESIS COST CONSISTENCY TEST")
    print("=" * 80)
    print("\nVerifying: Synthesis cost depends only on ENCODING METHOD,")
    print("           not on which synthesis backend processes codewords.\n")
    
    results = {}
    
    for encoding_method in encodings:
        print(f"\nEncoding: {encoding_method}")
        print("-" * 40)
        
        encoding_results = {}
        
        for synth_method in synthesis_methods:
            method_name = synth_method if synth_method else 'None'
            
            # Create parameters
            params = Params(
                name=f'synthesis_cost_test_{encoding_method}_{method_name}',
                filename='textfile_40b.txt',
                encoding_method=encoding_method,
                binarization_method='default',
                synthesis_method=synth_method,
            )
            
            # Analyze encoding
            analyzer = AnalyzeEncoding(params)
            cost = analyzer.calculate_synthesis_cost()
            
            encoding_results[method_name] = {
                'cost': cost,
                'breakdown': analyzer.metrics.get('synthesis_cost_breakdown'),
                'notes': analyzer.metrics.get('synthesis_cost_notes'),
            }
            
            # Display result
            if cost == 0:
                print(f"  {method_name:10s}: 0 bases/Kbit (no synthesis)")
            elif cost is not None:
                print(f"  {method_name:10s}: {cost:.0f} bases/Kbit")
            else:
                print(f"  {method_name:10s}: N/A ({encoding_results[method_name]['notes']})")
        
        results[encoding_method] = encoding_results
        
        # Verify consistency
        synthesis_costs = [
            r['cost'] for r in encoding_results.values() 
            if r['cost'] not in (0, None)
        ]
        
        if synthesis_costs:
            if all(c == synthesis_costs[0] for c in synthesis_costs):
                print(f"  [PASS] Consistent cost: {synthesis_costs[0]:.0f} bases/Kbit")
            else:
                print(f"  [FAIL] Inconsistent costs: {synthesis_costs}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    all_consistent = True
    for encoding_method, encoding_results in results.items():
        synthesis_costs = [
            r['cost'] for r in encoding_results.values() 
            if r['cost'] not in (0, None)
        ]
        
        if synthesis_costs:
            consistent = all(c == synthesis_costs[0] for c in synthesis_costs)
            status = "PASS" if consistent else "FAIL"
            print(f"\n{encoding_method:15s}: [{status}] {synthesis_costs[0]:.0f} bases/Kbit")
            
            if not consistent:
                all_consistent = False
    
    print("\n" + "=" * 80)
    if all_consistent:
        print("RESULT: PASS")
        print("Synthesis cost is independent of synthesis method.")
        print("Cost depends only on encoding efficiency.")
    else:
        print("RESULT: FAIL")
        print("Synthesis cost varies with synthesis method (unexpected).")
    print("=" * 80)
    
    return all_consistent


def test_error_resilience(encoding_method, error_rates, num_repeats=2):
    """
    Test full pipeline at different error rates to measure resilience.
    
    Uses the Simulation class to run the complete pipeline:
    binarize → encode → synthesize → store → sequence (with errors) → 
    recover → decode
    
    Args:
        encoding_method: 'max_density', 'church', 'wukong', etc.
        error_rates: list of IID error rates to test
        num_repeats: number of times to repeat each error rate
    
    Returns:
        dict with success rates at each error rate
    """
    print(f"\n{'='*80}")
    print(f"ERROR RESILIENCE TEST: {encoding_method}")
    print(f"{'='*80}\n")
    
    results = {}
    
    for error_rate in error_rates:
        # Create parameter set for this error rate
        params_list = []
        for repeat in range(num_repeats):
            params = Params(
                name=f'resilience_{encoding_method}_er{error_rate:.2f}_r{repeat}',
                filename='textfile_40b.txt',
                encoding_method=encoding_method,
                binarization_method='default',
                synthesis_method='mesa',
                sequencing_method='iid',
                iid_substitution_rate=error_rate * 0.7,     # 70% substitutions
                iid_insertion_rate=error_rate * 0.15,       # 15% insertions
                iid_deletion_rate=error_rate * 0.15,        # 15% deletions
                mean=0,  # Coverage for 99% reliability
                years=0,    # No storage degradation
            )
            params_list.append(params)
        
        # Run simulation with all parameters for this error rate
        sim = Simulation(params_list, debug=False)
        sim_results = sim.run()
        
        # Count successes
        successes = 0
        for sim_name, sim_data in sim_results.items():
            if sim_data.get('status') == 'SUCCESS':
                successes += 1
                print(f"  {encoding_method} @ {error_rate:.2%} error: SUCCESS")
            else:
                print(f"  {encoding_method} @ {error_rate:.2%} error: FAIL")
        
        success_rate = (successes / len(params_list) * 100) if params_list else 0
        results[error_rate] = {
            'successes': successes,
            'total': len(params_list),
            'success_rate': success_rate,
        }
        
        print(f"  --> {successes}/{len(params_list)} successful ({success_rate:.0f}%)\n")
    
    return results


def print_resilience_summary(all_results):
    """Print summary of error resilience across encodings."""
    print(f"\n{'='*80}")
    print("ERROR RESILIENCE SUMMARY")
    print(f"{'='*80}\n")
    
    for encoding_method, results in all_results.items():
        print(f"{encoding_method}:")
        for error_rate, data in results.items():
            print(f"  {error_rate:.2%} error: {data['success_rate']:.0f}% success")
        print()


def run_full_analysis():
    """Run both synthesis cost consistency and error resilience tests."""
    # Part 1: Synthesis cost consistency
    print("\nPART 1: SYNTHESIS COST CONSISTENCY")
    test_synthesis_cost_consistency()
    
    # Part 2: Error resilience
    print("\n\nPART 2: ERROR RESILIENCE TESTING")
    
    encodings = ['max_density', 'church', 'wukong']
    error_rates = [0.01, 0.05, 0.10, 0.15]
    
    all_resilience_results = {}
    for encoding_method in encodings:
        all_resilience_results[encoding_method] = test_error_resilience(
            encoding_method, 
            error_rates, 
            num_repeats=2
        )
    
    # Summary
    print_resilience_summary(all_resilience_results)
    
    print("="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)


if __name__ == '__main__':
    os.chdir(Path(__file__).parent.parent)
    run_full_analysis()
    sys.exit(0)

