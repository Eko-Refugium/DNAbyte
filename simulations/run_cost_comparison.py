"""
Example script: Compare encodings by cost using the simulation pipeline.

This script demonstrates how to use the cost comparison framework to:
1. Normalize all encodings to the same total DNA cost
2. Run the full simulation pipeline for each
3. Measure error resistance at equal cost
"""

import os
import json
import logging
from datetime import datetime
from typing import List, Dict, Any

from simulations.cost_comparison import (
    CostComparisonRunner, 
    CostCalculator,
    create_test_data
)
from simulations.simulation import Simulation
from dnabyte.data_classes.base import Data
from dnabyte.params import Params


class CostNormalizedComparison:
    """
    Runs simulations to compare encodings at equal cost.
    """
    
    def __init__(self, output_dir: str = './simulations/cost_normalized_results'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        
        handler = logging.FileHandler(
            os.path.join(output_dir, f'comparison_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
        )
        handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        self.logger.addHandler(handler)
    
    def run(self, 
            encodings: List[str] = None,
            test_file: str = './tests/testfiles/textfile_40b.txt',
            num_runs: int = 1,
            reference_encoding: str = 'goldman',
            with_errors: bool = False):
        """
        Run cost-normalized comparison using working parameters.
        
        Args:
            encodings: List of encoding names. If None, uses all 8 encodings.
            test_file: Test file to use
            num_runs: Number of simulation runs per encoding
            reference_encoding: Encoding to use as cost reference
            with_errors: If True, include sequencing errors
        """
        
        if encodings is None:
            encodings = ['goldman', 'church', 'gcplus', 'hedges', 
                        'max_density', 'no_homopolymer', 'wukong', 'yinyang']
        
        self.logger.info("=" * 80)
        self.logger.info("COST-NORMALIZED ENCODING COMPARISON")
        self.logger.info("=" * 80)
        
        # Working parameters for each encoding
        working_params = {
            'goldman': {
                'sequence_length': 200,
                'add_primer': False,
                'primer_length': 0,
            },
            'church': {
                'sequence_length': 200,
                'rs_num': 10,
                'add_primer': False,
                'primer_length': 0,
            },
            'gcplus': {
                'sequence_length': 200,
                'gcplus_k': 150,
                'add_primer': False,
                'primer_length': 0,
            },
            'hedges': {
                'sequence_length': 200,
                'hedges_coderate': 4,
                'add_primer': False,
                'primer_length': 0,
            },
            'max_density': {
                'codeword_length': 500,
                'add_primer': False,
                'primer_length': 0,
            },
            'no_homopolymer': {
                'codeword_length': 500,
                'max_homopolymer': 4,
                'add_primer': False,
                'primer_length': 0,
            },
            'wukong': {
                'sequence_length': 200,
                'rs_num': 10,
                'add_primer': False,
                'primer_length': 0,
            },
            'yinyang': {
                'sequence_length': 200,
                'max_content': 0.55,
                'add_primer': False,
                'primer_length': 0,
            },
        }
        
        # Load test data
        self.logger.info(f"Loading test data from {test_file}")
        test_data = Data(file_paths=[test_file])
        
        # Step 1: Run cost analysis
        self.logger.info("\n" + "=" * 80)
        self.logger.info("STEP 1: Cost Analysis")
        self.logger.info("=" * 80)
        
        cost_runner = CostComparisonRunner(
            output_dir=os.path.join(self.output_dir, 'cost_analysis')
        )
        
        # Build base_params from working params
        base_params = {
            'encoding_method': None,
            'binarization_method': 'default',
            'synthesis_method': 'nosynthpoly',
            'storage_conditions': None,
        }
        
        cost_results = cost_runner.run_cost_analysis(
            encodings=encodings,
            test_data=test_data,
            base_params=base_params,
            reference_encoding=reference_encoding,
            working_params=working_params
        )
        
        # Save cost analysis
        cost_summary = {
            enc: {
                'num_strands': metrics.num_strands,
                'avg_strand_length': metrics.avg_strand_length,
                'total_dna_length': metrics.total_dna_length,
                'parameters': metrics.parameters
            }
            for enc, metrics in cost_results.items()
        }
        
        with open(os.path.join(self.output_dir, 'cost_analysis.json'), 'w') as f:
            json.dump(cost_summary, f, indent=2)
        
        self.logger.info(f"Cost analysis saved to cost_analysis.json")
        
        # Step 2: Run full simulations for each encoding at normalized cost
        self.logger.info("\n" + "=" * 80)
        self.logger.info("STEP 2: Full Simulation Pipeline")
        self.logger.info("=" * 80)
        
        simulation_results = {}
        
        for enc in encodings:
            self.logger.info(f"\n--- Simulating {enc} ---")
            
            if enc not in cost_results:
                self.logger.warning(f"Skipping {enc}: cost analysis failed")
                continue
            
            # Create simulation parameters for this encoding
            params_list = []
            
            for run_id in range(num_runs):
                # Build Params object with working parameters
                params_dict = {
                    'name': f'cost_comp_{enc}_run{run_id + 1}',
                    'encoding_method': enc,
                    'binarization_method': 'default',
                    'add_primer': False,
                    'primer_length': 0,
                    'synthesis_method': 'nosynthpoly',
                    'storage_conditions': None,
                }
                
                # Add encoding-specific parameters
                if enc in working_params:
                    params_dict.update(working_params[enc])
                
                # Add error parameters if requested
                if with_errors:
                    params_dict.update({
                        'mean': 20,
                        'std_dev': 1,
                        'sequencing_method': 'kmere',
                        'kmer_k': 1,
                        'kmer_p_ins': 0.001,
                        'kmer_p_del': 0.001,
                        'kmer_p_sub': 0.002,
                        'kmer_seed': 42,
                    })
                
                params = Params(**params_dict)
                params_list.append(params)
            
            # Run simulation
            try:
                sim = Simulation(params_list, debug=False)
                results = sim.run(paralel=False)
                simulation_results[enc] = results
                
                self.logger.info(f"✓ Completed {num_runs} runs for {enc}")
                
                # Print summary
                success_count = sum(
                    1 for run_results in results.values() 
                    if isinstance(run_results, dict) and run_results.get('status') == 'SUCCESS'
                )
                self.logger.info(f"  Success rate: {success_count}/{num_runs}")
                
            except Exception as e:
                self.logger.error(f"✗ Simulation failed for {enc}: {str(e)}")
                import traceback
                self.logger.error(traceback.format_exc())
        
        # Step 3: Analyze and compare results
        self.logger.info("\n" + "=" * 80)
        self.logger.info("STEP 3: Comparison Summary")
        self.logger.info("=" * 80)
        
        self._print_comparison_summary(cost_results, simulation_results)
        
        # Save simulation results
        with open(os.path.join(self.output_dir, 'simulation_results.json'), 'w') as f:
            json.dump(simulation_results, f, indent=2, default=str)
        
        self.logger.info(f"\nResults saved to {self.output_dir}")
        
        return cost_results, simulation_results
    
    def _print_comparison_summary(self, cost_results: Dict, simulation_results: Dict):
        """Print a summary comparison of encodings."""
        
        self.logger.info("\nEncoding Costs:")
        self.logger.info("-" * 80)
        
        for enc in sorted(cost_results.keys()):
            metrics = cost_results[enc]
            self.logger.info(
                f"  {enc:20s}: "
                f"strands={metrics.num_strands:5d}, "
                f"avg_len={metrics.avg_strand_length:6.1f}, "
                f"total_dna={metrics.total_dna_length:7d}"
            )
        
        self.logger.info("\nError Resistance (from simulation):")
        self.logger.info("-" * 80)
        
        for enc in sorted(simulation_results.keys()):
            results = simulation_results[enc]
            
            success_count = sum(
                1 for run_results in results.values() 
                if run_results.get('status') == 'SUCCESS'
            )
            total_runs = len(results)
            success_rate = success_count / total_runs if total_runs > 0 else 0
            
            self.logger.info(
                f"  {enc:20s}: {success_rate*100:5.1f}% success "
                f"({success_count}/{total_runs} runs)"
            )


def main():
    """Run cost-normalized comparison for standard encodings."""
    
    encodings = ['goldman', 'church', 'gcplus', 'hedges', 'max_density', 'no_homopolymer', 'wukong', 'yinyang']
    
    comparison = CostNormalizedComparison()
    
    # Run with basic parameters (no errors)
    print("\n" + "=" * 80)
    print("RUNNING COST COMPARISON WITH WORKING PARAMETERS")
    print("=" * 80)
    
    comparison.run(
        encodings=encodings,
        test_file='./tests/testfiles/textfile_40b.txt',
        num_runs=1,
        reference_encoding='goldman',
        with_errors=False
    )


if __name__ == '__main__':
    main()
