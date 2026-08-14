#!/usr/bin/env python
"""
Quick-start script for cost-normalized encoding comparison.

Usage:
    python simulations/quickstart_cost_comparison.py

This script demonstrates the complete workflow:
1. Create test data
2. Analyze baseline costs for all 8 encodings
3. Compare error resistance at equal cost
"""

import os
import sys
import logging
from datetime import datetime

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Run quick comparison with working parameters."""
    
    print("\n" + "=" * 80)
    print("COST COMPARISON - QUICK START")
    print("=" * 80)
    
    from simulations.cost_comparison import CostComparisonRunner
    from dnabyte.data_classes.base import Data
    
    # Use larger test file for better gcplus divisibility (1950 bytes = 15600 bits)
    test_file = './tests/testfiles/Bohemian_Rhapsody_Lyrics.txt'
    
    if not os.path.exists(test_file):
        print(f"Error: Test file not found at {test_file}")
        print("Please run from the DNAbyte root directory")
        return
    
    # Create cost runner
    output_dir = './simulations/quickstart_results'
    runner = CostComparisonRunner(output_dir=output_dir)
    
    # Define encodings to test - all 8 encodings
    encodings = [
        'goldman', 'church', 'gcplus', 'hedges',
        'max_density', 'no_homopolymer', 'wukong', 'yinyang'
    ]
    
    # Working parameters for each encoding
    working_params = {
        'goldman': {'sequence_length': 200, 'add_primer': False, 'primer_length': 0},
        'church': {'sequence_length': 200, 'rs_num': 10, 'add_primer': False, 'primer_length': 0},
        'gcplus': {'sequence_length': 200, 'gcplus_k': 8, 'add_primer': False, 'primer_length': 0},  # 15600 bits divisible by 8
        'hedges': {'sequence_length': 200, 'hedges_coderate': 4, 'add_primer': False, 'primer_length': 0},
        'max_density': {'codeword_length': 200, 'add_primer': False, 'primer_length': 0},
        'no_homopolymer': {'codeword_length': 200, 'max_homopolymer': 4, 'add_primer': False, 'primer_length': 0},
        'wukong': {'sequence_length': 200, 'rs_num': 10, 'add_primer': False, 'primer_length': 0},
        'yinyang': {'sequence_length': 20, 'add_primer': False, 'primer_length': 0},  # Very small for performance
    }
    
    # Base parameters
    base_params = {
        'encoding_method': None,
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'storage_conditions': None,
    }
    
    # Load test data
    test_data = Data(file_paths=[test_file])
    
    print(f"\nTest file: {test_file}")
    print(f"Test data: {test_data.size * 8} bits ({test_data.size} bytes)")
    print(f"Encodings: {', '.join(encodings)}\n")
    
    # Run cost analysis
    print("Running cost analysis...")
    try:
        results = runner.run_cost_analysis(
            encodings=encodings,
            test_data=test_data,
            base_params=base_params,
            reference_encoding='goldman',
            working_params=working_params
        )
        
        print("\n" + "=" * 80)
        print("COST ANALYSIS RESULTS")
        print("=" * 80)
        
        for enc in sorted(results.keys()):
            metrics = results[enc]
            print(
                f"\n{enc}:\n"
                f"  Strands:        {metrics.num_strands}\n"
                f"  Avg length:     {metrics.avg_strand_length:.1f} bp\n"
                f"  Total DNA:      {metrics.total_dna_length} bp\n"
                f"  Bits/nucleotide: {metrics.bits_per_nucleotide:.4f}"
            )
        
        print(f"\nResults saved to: {output_dir}")
        
    except Exception as e:
        print(f"Error running cost analysis: {e}")
        import traceback
        traceback.print_exc()
        return


if __name__ == '__main__':
    main()
    results = None
    
    # Run advanced tuning example
    try:
        if results is not None:
            run_advanced_tuning_example()
    except Exception as e:
        logger.error(f"Advanced tuning failed: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
    
    logger.info("\n" + "=" * 80)
    logger.info("Quickstart complete!")
    logger.info("=" * 80)
    logger.info("\nNext steps:")
    logger.info("1. Review results in ./simulations/quickstart_results/")
    logger.info("2. Read COST_COMPARISON_GUIDE.md for detailed usage")
    logger.info("3. Modify run_cost_comparison.py for full simulation pipeline")
    logger.info("4. Use advanced_tuning.py for more sophisticated parameter optimization")


if __name__ == '__main__':
    main()
