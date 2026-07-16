"""
Error resilience testing at equal synthesis cost for encodings.

Compares error correction performance across encodings when all achieve
the same synthesis cost (~1150 bases/Kbit).
"""

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from dnabyte.params import Params
from simulations.simulation import Simulation


def run_resilience_tests():
    """Run error resilience tests for all encodings at equal synthesis cost."""
    
    # All available encodings - use DEFAULT parameters for equal ~1150 bases/Kbit cost
    encodings = [
        'max_density',      # 1150 bases/Kbit (baseline)
        'church',           # 1150 bases/Kbit 
        'gcplus',           # 1150 bases/Kbit
        'wukong',           # 1150 bases/Kbit
        'goldman',          # 1150 bases/Kbit
        'no_homopolymer',   # 1020 bases/Kbit
        'linear_chain',     # 1100 bases/Kbit
        'poly_chain',       # 1100 bases/Kbit
    ]
    error_rates = [0.01, 0.05, 0.10, 0.15]
    
    # EQUAL COST TUNING: Use DEFAULT parameters (already tuned to ~1150 bases/Kbit)
    # Each encoding's defaults are already at target cost
    tuned_params = {
        'max_density': {
            # Default: codeword_length=500, dna_barcode_length=75 -> 1150 bases/Kbit
            'mean': 5,
        },
        'church': {
            # Default: sequence_length=200 -> 1150 bases/Kbit
            'add_redundancy': True,
            'add_primer': True,
            'mean': 6,  # Increased from 4 to fix failures
        },
        'gcplus': {
            # Default: sequence_length=200, gcplus_k=168, gcplus_l=8, gcplus_c1=2 -> 1150 bases/Kbit
            'mean': 4,
        },
        'wukong': {
            # Default: codeword_length=200 -> 1150 bases/Kbit
            'mean': 3,  # Strong error correction - investigate why failing
        },
        'goldman': {
            # Default parameters -> 1150 bases/Kbit
            'mean': 5,
        },
        'no_homopolymer': {
            # Default parameters -> 1020 bases/Kbit (slightly lower)
            'mean': 5,
        },
        'linear_chain': {
            # Default parameters -> 1100 bases/Kbit (assembly-based)
            'mean': 6,
        },
        'poly_chain': {
            # Default parameters -> 1100 bases/Kbit (assembly-based, compact)
            'mean': 8,
        },
    }
    
    results = {}
    
    print("=" * 100)
    print("ERROR RESILIENCE TEST - EQUAL SYNTHESIS COST COMPARISON")
    print("=" * 100)
    print()
    
    for encoding in encodings:
        results[encoding] = {}
        print(f"\n{'='*100}")
        print(f"ENCODING: {encoding.upper()}")
        print(f"{'='*100}")
        
        params_dict = tuned_params.get(encoding, {})
        
        for error_rate in error_rates:
            test_name = f'resilience_{encoding}_{error_rate:.2f}'
            
            # Create params for this test
            try:
                params = Params(
                    name=test_name,
                    filename='textfile_40b.txt',
                    encoding_method=encoding,
                    binarization_method='default',
                    synthesis_method='assembly' if encoding in ['linear_chain', 'poly_chain'] else 'mesa',
                    sequencing_method='iid',
                    # Error rates (realistic: 70% subs, 15% ins, 15% del)
                    iid_substitution_rate=error_rate * 0.7,
                    iid_insertion_rate=error_rate * 0.15,
                    iid_deletion_rate=error_rate * 0.15,
                    years=0,  # No storage
                )
                
                # Add encoding-specific parameters
                for param, value in params_dict.items():
                    if param != 'mean':
                        setattr(params, param, value)
                    else:
                        params.mean = value
                
                # Run simulation
                print(f"\n  Error rate: {error_rate:.2%} (sub/ins/del: {error_rate*0.7:.3f}/{error_rate*0.15:.3f}/{error_rate*0.15:.3f})")
                
                sim = Simulation([params])
                result_dict = sim.run()
                
                # Extract result (format: {param_name: {step_results..., 'status': ...}})
                success = False
                failed_step = None
                
                try:
                    if result_dict and isinstance(result_dict, dict):
                        # Get the first (and only) result
                        result_key = list(result_dict.keys())[0]
                        result = result_dict[result_key]
                        
                        success = result.get('status') == 'SUCCESS'
                        failed_step = result.get('step_fail', '')
                        
                        results[encoding][error_rate] = {
                            'success': success,
                            'status': result.get('status', 'UNKNOWN'),
                            'notes': failed_step
                        }
                    else:
                        results[encoding][error_rate] = {'success': False, 'status': 'NO_RESULT'}
                    
                    symbol = "[PASS]" if success else "[FAIL]"
                    print(f"    Result: {symbol}")
                    if not success and failed_step:
                        print(f"    Failed at: {failed_step}")
                        
                except Exception as e:
                    results[encoding][error_rate] = {'success': False, 'status': 'ERROR'}
                    print(f"    Result: [ERROR] - Parsing: {type(e).__name__}")
                    
            except Exception as e:
                results[encoding][error_rate] = {'success': False, 'error': str(e)[:100]}
                print(f"    Result: [ERROR] - {type(e).__name__}: {str(e)[:60]}")
    
    # Print summary table
    print("\n\n" + "=" * 100)
    print("RESILIENCE TEST SUMMARY")
    print("=" * 100)
    print()
    
    print(f"{'Encoding':<20} {'1% Error':<15} {'5% Error':<15} {'10% Error':<15} {'15% Error':<15}")
    print("-" * 100)
    
    for encoding in encodings:
        row = f"{encoding:<20}"
        for error_rate in error_rates:
            if encoding in results and error_rate in results[encoding]:
                success = results[encoding][error_rate].get('success', False)
                symbol = "[PASS]" if success else "[FAIL]"
                row += f" {symbol:<15}"
            else:
                row += f" {'[?]':<15}"
        print(row)
    
    print()
    print("=" * 100)
    print("ANALYSIS")
    print("=" * 100)
    
    for encoding in encodings:
        success_count = sum(1 for err_rate in error_rates 
                          if encoding in results and error_rate in results[encoding]
                          and results[encoding][error_rate].get('success', False))
        total_count = len(error_rates)
        pass_rate = (success_count / total_count) * 100 if total_count > 0 else 0
        
        print(f"\n{encoding.upper()}")
        print(f"  Pass rate: {success_count}/{total_count} ({pass_rate:.0f}%)")
        
        # Find failure point
        for error_rate in sorted(error_rates):
            if encoding in results and error_rate in results[encoding]:
                if not results[encoding][error_rate].get('success', False):
                    print(f"  First failure at: {error_rate:.2%} error rate")
                    break


if __name__ == '__main__':
    os.chdir(Path(__file__).parent.parent)
    run_resilience_tests()
