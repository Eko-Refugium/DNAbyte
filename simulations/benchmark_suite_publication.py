"""
Comprehensive benchmark suite for DNAbyte publication
Runs all demonstration simulations and generates publication figures.

This script orchestrates the execution of multiple simulation studies that
demonstrate DNAbyte's capabilities for modular benchmarking of DNA storage systems.
"""
import os
import sys
from datetime import datetime
import logging

# Import the simulation runners
from simulations.sim_encoding_strategies_comparison import run_encoding_comparison
from simulations.sim_sequencing_error_impact import run_sequencing_error_impact
from simulations.sim_storage_encoding_interaction import run_storage_encoding_interaction
from simulations.sim_recovery_procedure_comparison import run_recovery_procedure_comparison


def setup_logging():
    """Set up logging for the benchmark suite."""
    log_dir = os.path.join('simulations', 'simlogs')
    os.makedirs(log_dir, exist_ok=True)
    
    log_filename = os.path.join(log_dir, f'benchmark_suite_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ]
    )
    
    return logging.getLogger(__name__)


def run_benchmark_suite():
    """
    Execute all benchmark simulations for the DNAbyte paper.
    
    Simulations demonstrate:
    1. Modular architecture - encoding strategies compared under identical conditions
    2. Configuration-driven simulation - parameter sweeps generating families of experiments
    3. Standardized benchmarking - reproducible comparisons across methodologies
    4. End-to-end evaluation - showing interactions between pipeline components
    """
    
    logger = setup_logging()
    
    logger.info("="*70)
    logger.info("DNAbyte Publication Benchmark Suite")
    logger.info("="*70)
    logger.info("Starting comprehensive benchmarking of DNA storage systems")
    
    # Track which simulations completed successfully
    completed = {}
    failed = {}
    
    # Simulation 1: Encoding strategy comparison
    logger.info("\n[1/4] Running: Encoding Strategy Comparison")
    logger.info("-" * 70)
    logger.info("Objective: Demonstrate modular architecture - compare encoding methods")
    logger.info("          with identical synthesis/storage/sequencing conditions")
    logger.info("Expected outcome: Success rates vary by encoding method")
    
    try:
        results1, agg1, rates1 = run_encoding_comparison()
        completed['encoding_comparison'] = {
            'results': results1,
            'aggregated': agg1,
            'success_rates': rates1
        }
        logger.info("[OK] Encoding comparison completed successfully")
    except Exception as e:
        logger.error(f"[FAIL] Encoding comparison failed: {str(e)}")
        failed['encoding_comparison'] = str(e)
    
    # Simulation 2: Sequencing error impact
    logger.info("\n[2/4] Running: Sequencing Error Impact Analysis")
    logger.info("-" * 70)
    logger.info("Objective: Show how sequencing technology affects encoding performance")
    logger.info("          across realistic error rate ranges (0.1% to 5%)")
    logger.info("Expected outcome: Different encodings show different error sensitivity")
    
    try:
        results2, agg2, rates2 = run_sequencing_error_impact()
        completed['sequencing_error_impact'] = {
            'results': results2,
            'aggregated': agg2,
            'success_data': rates2
        }
        logger.info("[OK] Sequencing error impact completed successfully")
    except Exception as e:
        logger.error(f"[FAIL] Sequencing error impact failed: {str(e)}")
        failed['sequencing_error_impact'] = str(e)
    
    # Simulation 3: Storage-encoding interaction
    logger.info("\n[3/4] Running: Storage Duration and Encoding Interaction")
    logger.info("-" * 70)
    logger.info("Objective: Demonstrate end-to-end interactions between storage")
    logger.info("          conditions and encoding strategy choice")
    logger.info("          across 1-10,000 year timescales")
    logger.info("Expected outcome: Encoding performance degrades differently with time")
    
    try:
        results3, agg3, rates3 = run_storage_encoding_interaction()
        completed['storage_encoding_interaction'] = {
            'results': results3,
            'aggregated': agg3,
            'success_data': rates3
        }
        logger.info("[OK] Storage-encoding interaction completed successfully")
    except Exception as e:
        logger.error(f"[FAIL] Storage-encoding interaction failed: {str(e)}")
        failed['storage_encoding_interaction'] = str(e)
    
    # Simulation 4: Recovery procedure comparison
    logger.info("\n[4/4] Running: Recovery Procedure Comparison")
    logger.info("-" * 70)
    logger.info("Objective: Compare different clustering+recovery method combinations")
    logger.info("          showing impact of recovery strategy on performance")
    logger.info("Expected outcome: Recovery method choice affects robustness to errors")
    
    try:
        results4, agg4, rates4 = run_recovery_procedure_comparison()
        completed['recovery_comparison'] = {
            'results': results4,
            'aggregated': agg4,
            'success_data': rates4
        }
        logger.info("[OK] Recovery procedure comparison completed successfully")
    except Exception as e:
        logger.error(f"[FAIL] Recovery procedure comparison failed: {str(e)}")
        failed['recovery_comparison'] = str(e)
    
    # Summary
    logger.info("\n" + "="*70)
    logger.info("BENCHMARK SUITE SUMMARY")
    logger.info("="*70)
    logger.info(f"Completed: {len(completed)}/4 simulations")
    for sim_name in completed.keys():
        logger.info(f"  [OK] {sim_name}")
    
    if failed:
        logger.info(f"Failed: {len(failed)}/4 simulations")
        for sim_name, error in failed.items():
            logger.error(f"  [FAIL] {sim_name}: {error}")
    
    logger.info("\n" + "-"*70)
    logger.info("All benchmark figures have been generated and saved to:")
    logger.info(f"  {os.path.join('simulations', 'simlogs')}")
    logger.info("\nFigures generated:")
    logger.info("  1. encoding_comparison_*.png")
    logger.info("     Shows comparative performance of different encoding strategies")
    logger.info("  2. sequencing_error_impact_*.png")
    logger.info("     Shows how sequencing errors differentially affect encodings")
    logger.info("  3. storage_encoding_interaction_*.png")
    logger.info("     Shows storage duration effects on encoding performance")
    logger.info("  4. recovery_comparison_*.png")
    logger.info("     Shows impact of recovery procedure choice")
    logger.info("\nUse these figures in:")
    logger.info("  - Paper Figures (Results section)")
    logger.info("  - Supplementary material for benchmark details")
    logger.info("  - Documentation of framework capabilities")
    
    logger.info("\n" + "="*70)
    logger.info("Benchmark suite execution completed")
    logger.info("="*70)
    
    return completed, failed


if __name__ == '__main__':
    completed, failed = run_benchmark_suite()
    
    # Exit with appropriate code
    sys.exit(0 if len(failed) == 0 else 1)
