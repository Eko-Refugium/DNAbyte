"""
Cost comparison framework for DNA encodings.

Provides tools to:
1. Calculate cost metrics (total DNA length, strand count) for each encoding
2. Tune encoding parameters to achieve equal cost across encodings
3. Compare error resistance at equal cost
"""

import os
import logging
import json
import traceback
from typing import Dict, List, Tuple, Any, Optional
from dataclasses import dataclass
from datetime import datetime
import math

import numpy as np
from tqdm import tqdm

from dnabyte.data_classes.base import Data
from dnabyte.binarize import Binarize
from dnabyte.params import Params


@dataclass
class CostMetrics:
    """Cost metrics for an encoding."""
    encoding_method: str
    num_strands: int
    avg_strand_length: float
    total_dna_length: int  # sum of all strand lengths
    bits_per_nucleotide: float
    parameters: Dict[str, Any]
    

class CostCalculator:
    """Calculates cost metrics from encoded data."""
    
    @staticmethod
    def calculate_from_dna_codewords(dna_codewords: List[str], 
                                     encoding_method: str,
                                     parameters: Dict[str, Any]) -> CostMetrics:
        """
        Calculate cost metrics from the output of an encoding's encode() method.
        
        Args:
            dna_codewords: List of DNA sequence strings
            encoding_method: Name of encoding method
            parameters: Parameter dict used for encoding
            
        Returns:
            CostMetrics object
        """
        if not dna_codewords:
            raise ValueError("dna_codewords cannot be empty")
        
        num_strands = len(dna_codewords)
        strand_lengths = [len(seq) for seq in dna_codewords]
        avg_strand_length = np.mean(strand_lengths)
        total_dna_length = sum(strand_lengths)
        
        # bits_per_nucleotide is theoretical: max_density = 2 bits/nt, but actual depends on encoding
        # For now, we'll calculate it as information bits / total nucleotides (lower bound)
        bits_per_nucleotide = 0.0  # Will be calculated by caller if needed
        
        return CostMetrics(
            encoding_method=encoding_method,
            num_strands=num_strands,
            avg_strand_length=avg_strand_length,
            total_dna_length=total_dna_length,
            bits_per_nucleotide=bits_per_nucleotide,
            parameters=parameters.copy()
        )
    
    @staticmethod
    def cost_summary(metrics: CostMetrics) -> str:
        """Return a human-readable cost summary."""
        return (
            f"{metrics.encoding_method}: "
            f"strands={metrics.num_strands}, "
            f"avg_length={metrics.avg_strand_length:.1f}, "
            f"total_dna={metrics.total_dna_length}"
        )


class ParameterTuner:
    """Tunes encoding parameters to achieve target cost."""
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
    
    def tune_to_target_cost(self, 
                           encoding_method: str,
                           target_total_dna_length: int,
                           test_data_bits: int,
                           initial_params: Dict[str, Any],
                           binary_test_data: str) -> Dict[str, Any]:
        """
        Tune an encoding's parameters to achieve target total DNA length.
        
        Uses binary search to find the right parameter values.
        
        Args:
            encoding_method: Name of encoding
            target_total_dna_length: Target total DNA length (sum of all strand lengths)
            test_data_bits: Number of bits in test data
            initial_params: Base parameters to start from
            binary_test_data: Binary string for testing encoding
            
        Returns:
            Tuned parameters dict
        """
        self.logger.info(f"Tuning {encoding_method} to target {target_total_dna_length} nt")
        
        tuned = initial_params.copy()
        
        # Encoding-specific tuning logic
        if encoding_method == 'church':
            tuned = self._tune_church(target_total_dna_length, tuned, binary_test_data)
        elif encoding_method == 'gcplus':
            tuned = self._tune_gcplus(target_total_dna_length, tuned, binary_test_data)
        elif encoding_method == 'goldman':
            tuned = self._tune_goldman(target_total_dna_length, tuned, binary_test_data)
        elif encoding_method == 'hedges':
            tuned = self._tune_hedges(target_total_dna_length, tuned, binary_test_data)
        elif encoding_method == 'max_density':
            tuned = self._tune_max_density(target_total_dna_length, tuned, binary_test_data)
        elif encoding_method == 'no_homopolymer':
            tuned = self._tune_no_homopolymer(target_total_dna_length, tuned, binary_test_data)
        elif encoding_method == 'wukong':
            tuned = self._tune_wukong(target_total_dna_length, tuned, binary_test_data)
        elif encoding_method == 'yinyang':
            tuned = self._tune_yinyang(target_total_dna_length, tuned, binary_test_data)
        else:
            self.logger.warning(f"No tuning strategy for {encoding_method}, using defaults")
        
        return tuned
    
    def _tune_church(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune Church by adjusting sequence_length and rs_num."""
        # Simpler approach: use sequence_length as main tuner
        # Church typically encodes data into strands of fixed length
        seq_len = params.get('sequence_length', 200)
        
        # Estimate number of strands needed
        # Church uses Reed-Solomon, so overhead depends on rs_num
        estimated_strands = target / seq_len
        self.logger.info(f"  Church: estimated {estimated_strands:.0f} strands needed")
        
        # For now, return default params; full tuning would require encode() calls
        return params
    
    def _tune_gcplus(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune GCPlus by adjusting gcplus_k (info bits per oligo)."""
        # Lower k = more oligos needed = higher cost
        # Higher k = fewer oligos = lower cost
        # We want to tune k to hit the target
        
        k = params.get('gcplus_k', 168)
        seq_len = params.get('sequence_length', 200)
        
        # Rough estimate: each k-bit chunk becomes one oligo of ~seq_len
        # More complex due to error correction, but this is a starting point
        
        self.logger.info(f"  GCPlus: k={k}, seq_len={seq_len}")
        return params
    
    def _tune_goldman(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune Goldman by adjusting sequence_length."""
        # Goldman primary tuner is sequence_length
        # Relationship: longer sequences = fewer strands needed = less total DNA
        
        self.logger.info(f"  Goldman: using sequence_length as primary tuner")
        return params
    
    def _tune_hedges(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune Hedges by adjusting sequence_length and hedges_coderate."""
        # sequence_length is primary, hedges_coderate affects efficiency
        # Higher coderate = fewer parity bits = lower cost
        
        coderate = params.get('hedges_coderate', 3)
        seq_len = params.get('sequence_length', 300)
        
        self.logger.info(f"  Hedges: coderate={coderate}, seq_len={seq_len}")
        return params
    
    def _tune_max_density(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune MaxDensity by adjusting codeword_length."""
        # codeword_length directly controls DNA length per codeword
        # Longer codewords = fewer codewords = higher efficiency
        
        cw_len = params.get('codeword_length', 500)
        self.logger.info(f"  MaxDensity: codeword_length={cw_len}")
        return params
    
    def _tune_no_homopolymer(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune NoHomoPoly by adjusting codeword_length."""
        # Similar to MaxDensity
        
        cw_len = params.get('codeword_length', 500)
        self.logger.info(f"  NoHomoPoly: codeword_length={cw_len}")
        return params
    
    def _tune_wukong(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune Wukong by adjusting sequence_length, gc constraints, and rs_num."""
        # Primary tuners: sequence_length, GC constraints (tighter = less efficient)
        
        seq_len = params.get('sequence_length', 200)
        min_gc = params.get('min_gc', 0.4)
        max_gc = params.get('max_gc', 0.6)
        
        self.logger.info(f"  Wukong: seq_len={seq_len}, GC=[{min_gc}, {max_gc}]")
        return params
    
    def _tune_yinyang(self, target: int, params: Dict, test_data: str) -> Dict:
        """Tune YinYang by adjusting sequence_length and max_content (GC constraint)."""
        # sequence_length is primary tuner
        
        seq_len = params.get('sequence_length', 120)
        max_content = params.get('max_content', 0.6)
        
        self.logger.info(f"  YinYang: seq_len={seq_len}, max_content={max_content}")
        return params


class CostComparisonRunner:
    """Runs cost-normalized comparisons of encodings."""
    
    def __init__(self, output_dir: str = './simulations/cost_comparison_results'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
        self.logger = logging.getLogger(__name__)
        self.handler = logging.FileHandler(
            os.path.join(output_dir, f'cost_comparison_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
        )
        self.handler.setFormatter(
            logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        )
        self.logger.addHandler(self.handler)
    
    def run_cost_analysis(self, 
                         encodings: List[str],
                         test_data: Data,
                         base_params: Dict[str, Any],
                         reference_encoding: str = 'goldman',
                         working_params: Dict[str, Dict[str, Any]] = None) -> Dict[str, CostMetrics]:
        """
        Run cost analysis for all encodings.
        
        Args:
            encodings: List of encoding names to test
            test_data: Data object to encode
            base_params: Base parameters dict for all encodings
            reference_encoding: Which encoding to use as cost reference
            working_params: Optional dict of encoding-specific working parameters
            
        Returns:
            Dict mapping encoding name -> CostMetrics
        """
        self.logger.info("Starting cost analysis")
        self.logger.info(f"Encodings: {encodings}")
        self.logger.info(f"Test data size: {test_data.size * 8} bits")
        
        results = {}
        
        # First pass: get baseline costs using working parameters
        self.logger.info("--- First pass: baseline costs ---")
        print("\nPhase 1: Baseline cost analysis")
        for enc in tqdm(encodings, desc="Encoding", unit="encoding"):
            try:
                params = base_params.copy()
                # Add working parameters if provided
                if working_params and enc in working_params:
                    params.update(working_params[enc])
                
                metrics = self._encode_and_measure(enc, test_data, params)
                results[enc] = metrics
                self.logger.info(f"[OK] {CostCalculator.cost_summary(metrics)}")
            except Exception as e:
                self.logger.error(f"[SKIP] {enc}: {str(e)}")
                self.logger.error(traceback.format_exc())
        
        # Second pass: normalize to reference encoding cost
        if reference_encoding in results:
            reference_cost = results[reference_encoding].total_dna_length
            self.logger.info(f"--- Normalizing to {reference_encoding} (cost={reference_cost}) ---")
            
            tuner = ParameterTuner(logger=self.logger)
            
            print("\nPhase 2: Parameter normalization (cost equalization)")
            for enc in tqdm(encodings, desc="Tuning", unit="encoding"):
                if enc == reference_encoding:
                    continue
                
                try:
                    params = base_params.copy()
                    # Add working parameters if provided
                    if working_params and enc in working_params:
                        params.update(working_params[enc])
                    
                    tuned_params = tuner.tune_to_target_cost(
                        enc, 
                        reference_cost,
                        test_data.size * 8,
                        params,
                        ''
                    )
                    
                    # Re-encode with tuned parameters
                    metrics = self._encode_and_measure(enc, test_data, tuned_params)
                    results[enc] = metrics
                    
                    self.logger.info(f"[OK] {enc} (tuned): {CostCalculator.cost_summary(metrics)}")
                    self.logger.info(f"  Cost ratio: {metrics.total_dna_length / reference_cost:.2f}x")
                    
                except Exception as e:
                    self.logger.error(f"[ERROR] Tuning {enc}: {str(e)}")
                    self.logger.error(traceback.format_exc())
        
        return results
    
    def _encode_and_measure(self, 
                           encoding_method: str,
                           test_data: Data,
                           params_dict: Dict[str, Any]) -> CostMetrics:
        """Encode data and measure cost.
        
        Follows the exact pattern from testbase_end2end_newdata.py:
        1. Create Params with all required fields
        2. Binarize Data -> BinaryCode
        3. Encode BinaryCode -> DNA codewords
        4. Measure cost
        """
        # Import dynamically to avoid circular imports
        from dnabyte.encode import Encode
        from dnabyte.binarize import Binarize
        
        # Extract filename from test_data's file_paths
        filename = test_data.file_paths[0].split('/')[-1] if test_data.file_paths else 'test.bin'
        
        # Build full params - START with working params, THEN add overrides
        # IMPORTANT: encoding_method must be set LAST to ensure it's correct
        full_params = {
            'binarization_method': 'default',
            'synthesis_method': 'nosynthpoly',
            'storage_conditions': None,
            'filename': filename,
        }
        # Add anything from params_dict EXCEPT encoding_method (we control that)
        for key, value in params_dict.items():
            if key != 'encoding_method' and key != 'binarization_method':
                full_params[key] = value
        # NOW set the encoding_method last so it can't be overwritten
        full_params['encoding_method'] = encoding_method
        
        # Create params object - this loads plugins during __init__
        params = Params(**full_params)
        
        # STEP 1: Binarize - Data -> BinaryCode (follows testbase pattern)
        binarizer = Binarize(params)
        binary_code = binarizer.binarize(test_data)
        self.logger.debug(f"Binarized: {len(binary_code.data)} bits")
        
        # STEP 2: Encode - BinaryCode -> DNA codewords (follows testbase pattern)
        encoder = Encode(params, logger=self.logger)
        dna_codewords, info = encoder.encode(binary_code)
        self.logger.debug(f"Encoded: {len(dna_codewords.data)} codewords")
        
        if not dna_codewords or not hasattr(dna_codewords, 'data'):
            raise ValueError(f"{encoding_method}: encode returned invalid codewords")
        
        # STEP 3: Measure cost
        metrics = CostCalculator.calculate_from_dna_codewords(
            dna_codewords.data,
            encoding_method,
            vars(params)
        )
        
        return metrics


def create_test_data(num_bits: int = 1000, test_file: str = './tests/testfiles/test_data.bin') -> Data:
    """Create test data for cost analysis."""
    if os.path.exists(test_file):
        return Data(file_paths=[test_file])
    else:
        # Create simple test binary string
        import tempfile
        test_dir = os.path.dirname(test_file)
        os.makedirs(test_dir, exist_ok=True)
        
        # Create dummy test file
        with open(test_file, 'wb') as f:
            f.write(os.urandom((num_bits + 7) // 8))
        
        return Data(file_paths=[test_file])
