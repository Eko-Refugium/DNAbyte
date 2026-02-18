"""
Standalone benchmark for assembly with parallelization.
This runs outside unittest context, allowing threading/multiprocessing.
"""
import time
import sys
import os

# Add DNAbyte to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dnabyte.binarize import binarize
from dnabyte.encode import encode
from dnabyte.store import synthesize, store, sequence
from dnabyte.params import Params

def benchmark_assembly(mean=500, years=100, encoding_method='poly_chain'):
    """Run assembly benchmark with large data."""
    
    print(f"\n{'='*60}")
    print(f"Benchmarking Assembly with mean={mean}, years={years}")
    print(f"{'='*60}\n")
    
    # Step 1: Binarization
    start_time = time.time()
    binarization_method = 'default'
    binarize_object = binarize(
        '../tests/testfiles/textfile_40b.txt',
        binarization_method=binarization_method,
        out_dir='tests/testdecode/'
    )
    print(f"✓ Binarization: {time.time() - start_time:.3f}s")
    
    # Step 2: Encoding
    start_time = time.time()
    library_file = 'polymeraselibinparts.txt'
    encode_object = encode(binarize_object, encoding_method=encoding_method, library=library_file)
    print(f"✓ Encoding: {time.time() - start_time:.3f}s")
    
    # Step 3: Assembly (THIS IS WHAT WE'RE BENCHMARKING)
    print(f"\n🔧 Starting Assembly with mean={mean}...")
    synthesis_start = time.time()
    
    params = Params(library_file=library_file, synthesis_method='assembly',
                   theory='no', mean=mean, std_dev=10, 
                   hybridisation_steps=10000)
    
    synthesized_object = synthesize(encode_object, params=params, use_parallel=True)
    
    assembly_time = time.time() - synthesis_start
    print(f"✓ Assembly completed: {assembly_time:.2f}s")
    
    # Step 4: Storage
    start_time = time.time()
    storage_params = Params(storage_method='biogene', temperature=37.0, time=years)
    stored_object = store(synthesized_object, params=storage_params)
    print(f"✓ Storage simulation: {time.time() - start_time:.3f}s")
    
    # Step 5: Sequencing
    start_time = time.time()
    sequencing_params = Params(sequencing_method='iid', iid_error_rate=0.05)
    sequenced_object = sequence(stored_object, params=sequencing_params)
    print(f"✓ Sequencing simulation: {time.time() - start_time:.3f}s")
    
    print(f"\n{'='*60}")
    print(f"ASSEMBLY TIME: {assembly_time:.2f}s")
    print(f"{'='*60}\n")
    
    return assembly_time


if __name__ == '__main__':
    # Run benchmarks with different parameters
    import multiprocessing as mp
    mp.set_start_method('spawn', force=True)
    
    benchmarks = [
        (20, 1, 'poly_chain'),
        (100, 10, 'poly_chain'),
        (500, 100, 'poly_chain'),
    ]
    
    results = []
    for mean, years, method in benchmarks:
        try:
            assembly_time = benchmark_assembly(mean, years, method)
            results.append((mean, years, assembly_time))
        except Exception as e:
            print(f"❌ Benchmark failed for mean={mean}: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary
    print(f"\n{'='*60}")
    print("BENCHMARK SUMMARY")
    print(f"{'='*60}")
    for mean, years, assembly_time in results:
        print(f"mean={mean:3d}, years={years:3d} → Assembly: {assembly_time:6.2f}s")
    print(f"{'='*60}\n")
