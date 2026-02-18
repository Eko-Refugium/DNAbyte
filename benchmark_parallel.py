"""
Benchmark script that runs OUTSIDE unittest to test parallel assembly with large parameters.
This demonstrates that parallelization works when not constrained by unittest's daemon context.
"""
import time
import sys
from scipy.constants import Avogadro

from dnabyte.data_classes.base import Data
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.params import Params

def benchmark_parallel_assembly():
    print("\n" + "="*70)
    print("PARALLEL ASSEMBLY BENCHMARK (Outside unittest)")
    print("="*70 + "\n")
    
    # Test with increasing mean values
    test_cases = [
        (50, "Small"),
        (100, "Medium"),
        (500, "Large - should be fast with parallelization!")
    ]
    
    for mean_val, description in test_cases:
        print(f"\n{'='*70}")
        print(f"Testing with mean={mean_val} ({description})")
        print('='*70)
        
        # Create params with current mean value
        params = Params(
            filename='textfile_40b.txt',
            encoding_method='poly_chain',
            library_name='polymeraselibinparts.txt',
            binarization_method='default',
            synthesis_method='assembly',
            theory='no',
            mean=mean_val,
            std_dev=0,
            hybridisation_steps=10000
        )
        
        # Load and binarize data
        print("Loading and binarizing data...")
        data = Data(file_paths=['tests/testfiles/textfile_40b.txt'])
        bin_obj = Binarize(params)
        binary_code = bin_obj.binarize(data)
        
        print("Encoding...")
        enc = Encode(params)
        data_enc, info = enc.encode(binary_code)
        
        print(f"Starting assembly with mean={mean_val}...")
        start = time.time()
        syn = SimulateSynthesis(params)
        data_syn, _ = syn.simulate(data_enc)
        elapsed = time.time() - start
        
        print(f"\n✓ ASSEMBLY COMPLETED in {elapsed:.2f}s")
        
        if mean_val == 500 and elapsed < 30:
            print("  ⚡ PARALLELIZATION IS WORKING! (much faster than sequential)")
        elif mean_val == 500 and elapsed < 60:
            print("  ⚡ Good speed - threading appears to be working")
        elif mean_val == 500:
            print("  ⚠ Still slow - may be running sequentially")
    
    print(f"\n{'='*70}")
    print("BENCHMARK COMPLETE")
    print("="*70 + "\n")

if __name__ == '__main__':
    # This script runs in normal Python context (not unittest daemon)
    # so parallelization should work
    print("\nNote: This runs outside unittest, so threading/multiprocessing works properly.")
    print("Compare this to unittest results to see the difference!\n")
    
    try:
        benchmark_parallel_assembly()
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
