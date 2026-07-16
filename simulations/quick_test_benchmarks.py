"""
Quick test version of benchmark suite - uses minimal repeats for fast verification
Checks that simulations run without errors before running full version
"""
import os
import logging
from datetime import datetime

from dnabyte.params import Params
from simulations.simulation import Simulation


def quick_test_encoding_comparison():
    """Quick test with 1 repeat per encoding (not 5)"""
    print("\nQUICK TEST: Encoding Strategy Comparison")
    print("-" * 50)
    
    encoding_methods = ['max_density']  # Just one for speed
    repeats = 1  # Just one repeat
    
    params_list = []
    for encoding_method in encoding_methods:
        for repeat in range(repeats):
            params = Params(
                name=f'quick_test_{encoding_method}_{repeat+1}',
                filename='textfile_40b.txt',
                encoding_method=encoding_method,
                binarization_method='default',
                synthesis_method=None,
                storage_conditions='biogene',
                sequencing_method=None,
                years=10,
                codeword_length=500,
                dna_barcode_length=75,
            )
            params_list.append(params)
    
    print(f"Running {len(params_list)} quick test simulations...")
    sim = Simulation(params_list)
    results = sim.run(paralel=False)
    
    # Check results
    success = sum(1 for r in results.values() if r.get('status') == 'SUCCESS')
    total = len(results)
    
    print(f"Results: {success}/{total} successful")
    if success == total:
        print("[OK] Quick test passed!")
        return True
    else:
        print("[FAIL] Some simulations failed")
        return False


def quick_test_all():
    """Run all quick tests"""
    print("\n" + "="*70)
    print("RUNNING QUICK VERIFICATION TESTS")
    print("="*70)
    
    all_passed = True
    
    try:
        if quick_test_encoding_comparison():
            print("\n[OK] Encoding comparison works")
        else:
            all_passed = False
            print("\n[FAIL] Encoding comparison failed")
    except Exception as e:
        print(f"\n[ERROR] Encoding comparison: {e}")
        all_passed = False
    
    print("\n" + "="*70)
    if all_passed:
        print("All quick tests PASSED - ready to run full benchmark suite")
        print("\nRun full suite with:")
        print("  python simulations/benchmark_suite_publication.py")
    else:
        print("Some quick tests FAILED - fix errors before running full suite")
    print("="*70)
    
    return all_passed


if __name__ == '__main__':
    success = quick_test_all()
    exit(0 if success else 1)
