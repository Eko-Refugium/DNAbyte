"""
Advanced parameter tuning with actual encoding calls.

This module implements sophisticated parameter tuning that uses binary search
and actual encoding to find the right parameters for cost-normalized comparisons.
"""

import logging
from typing import Dict, Any, Callable, Tuple, Optional
from functools import lru_cache

from dnabyte.params import Params
from dnabyte.encode import Encode
from dnabyte.data_classes.base import Data
from simulations.cost_comparison import CostCalculator, CostMetrics


class AdvancedParameterTuner:
    """
    Advanced parameter tuning using binary search and actual encoding calls.
    """
    
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self._encode_cache = {}  # Cache encoding results
    
    def tune_encoding_to_cost(self,
                             encoding_method: str,
                             target_cost: int,
                             test_data: Data,
                             base_params: Dict[str, Any],
                             tuning_params: Dict[str, Tuple[float, float]],
                             tolerance: int = 100,
                             max_iterations: int = 20) -> Dict[str, Any]:
        """
        Tune encoding parameters to achieve target cost using binary search.
        
        Args:
            encoding_method: Name of encoding
            target_cost: Target total DNA length in nucleotides
            test_data: Test data to encode
            base_params: Base parameters
            tuning_params: Dict of param_name -> (min_val, max_val) ranges for tuning
            tolerance: Acceptable cost deviation from target
            max_iterations: Max binary search iterations per parameter
            
        Returns:
            Tuned parameters dict
        """
        self.logger.info(f"Tuning {encoding_method} to cost {target_cost} nt (±{tolerance})")
        
        tuned = base_params.copy()
        
        # Tune each parameter sequentially
        for param_name, (min_val, max_val) in tuning_params.items():
            self.logger.info(f"  Tuning {param_name} in range [{min_val}, {max_val}]")
            
            tuned[param_name] = self._binary_search_parameter(
                encoding_method=encoding_method,
                param_name=param_name,
                min_val=min_val,
                max_val=max_val,
                target_cost=target_cost,
                test_data=test_data,
                base_params=tuned,
                tolerance=tolerance,
                max_iterations=max_iterations
            )
            
            current_cost = self._encode_and_measure_cost(
                encoding_method, test_data, tuned
            )
            self.logger.info(
                f"    → {param_name} = {tuned[param_name]}, cost = {current_cost}"
            )
        
        return tuned
    
    def _binary_search_parameter(self,
                                encoding_method: str,
                                param_name: str,
                                min_val: float,
                                max_val: float,
                                target_cost: int,
                                test_data: Data,
                                base_params: Dict[str, Any],
                                tolerance: int,
                                max_iterations: int) -> float:
        """
        Binary search to find optimal parameter value.
        
        Strategy: vary the parameter and measure resulting cost.
        """
        
        # Determine parameter range and step
        if isinstance(min_val, int) and isinstance(max_val, int):
            # Integer parameter
            search_fn = self._binary_search_int
        else:
            # Float parameter
            search_fn = self._binary_search_float
        
        best_value = search_fn(
            encoding_method=encoding_method,
            param_name=param_name,
            min_val=min_val,
            max_val=max_val,
            target_cost=target_cost,
            test_data=test_data,
            base_params=base_params,
            tolerance=tolerance,
            max_iterations=max_iterations
        )
        
        return best_value
    
    def _binary_search_int(self,
                          encoding_method: str,
                          param_name: str,
                          min_val: int,
                          max_val: int,
                          target_cost: int,
                          test_data: Data,
                          base_params: Dict[str, Any],
                          tolerance: int,
                          max_iterations: int) -> int:
        """Binary search for integer parameter."""
        
        lo, hi = int(min_val), int(max_val)
        best_value = lo
        best_cost = float('inf')
        
        for iteration in range(max_iterations):
            mid = (lo + hi) // 2
            
            # Test this value
            test_params = base_params.copy()
            test_params[param_name] = mid
            
            try:
                cost = self._encode_and_measure_cost(
                    encoding_method, test_data, test_params
                )
            except Exception as e:
                self.logger.debug(f"Failed to encode with {param_name}={mid}: {str(e)}")
                # Skip invalid parameter values
                if mid < target_cost:
                    lo = mid + 1
                else:
                    hi = mid - 1
                continue
            
            # Track best solution
            cost_error = abs(cost - target_cost)
            if cost_error < abs(best_cost - target_cost):
                best_value = mid
                best_cost = cost
            
            # Check convergence
            if cost_error <= tolerance:
                self.logger.debug(f"Converged at {param_name}={mid} (cost={cost})")
                return mid
            
            # Adjust search range
            if cost > target_cost:
                # Cost is too high, need to reduce it
                # Depends on parameter direction (some params increase cost when higher, some decrease)
                # For now, try both directions and pick the better one
                hi = mid - 1
            else:
                lo = mid + 1
            
            if lo > hi:
                break
        
        return best_value
    
    def _binary_search_float(self,
                            encoding_method: str,
                            param_name: str,
                            min_val: float,
                            max_val: float,
                            target_cost: int,
                            test_data: Data,
                            base_params: Dict[str, Any],
                            tolerance: int,
                            max_iterations: int) -> float:
        """Binary search for float parameter."""
        
        lo, hi = float(min_val), float(max_val)
        best_value = lo
        best_cost = float('inf')
        
        for iteration in range(max_iterations):
            mid = (lo + hi) / 2.0
            
            # Test this value
            test_params = base_params.copy()
            test_params[param_name] = mid
            
            try:
                cost = self._encode_and_measure_cost(
                    encoding_method, test_data, test_params
                )
            except Exception as e:
                self.logger.debug(f"Failed to encode with {param_name}={mid}: {str(e)}")
                # Try different direction
                if mid < 0.5 * (lo + hi):
                    lo = mid
                else:
                    hi = mid
                continue
            
            # Track best solution
            cost_error = abs(cost - target_cost)
            if cost_error < abs(best_cost - target_cost):
                best_value = mid
                best_cost = cost
            
            # Check convergence
            if cost_error <= tolerance:
                self.logger.debug(f"Converged at {param_name}={mid:.4f} (cost={cost})")
                return mid
            
            # Adjust search range
            if cost > target_cost:
                hi = mid
            else:
                lo = mid
            
            if hi - lo < 1e-6:
                break
        
        return best_value
    
    def _encode_and_measure_cost(self,
                                encoding_method: str,
                                test_data: Data,
                                params_dict: Dict[str, Any]) -> int:
        """
        Encode data and return total DNA cost.
        
        Returns:
            Total DNA length in nucleotides
        """
        
        # Check cache first
        cache_key = self._make_cache_key(encoding_method, params_dict)
        if cache_key in self._encode_cache:
            return self._encode_cache[cache_key]
        
        # Create params object
        params = Params(encoding_method=encoding_method)
        for key, value in params_dict.items():
            if key != 'encoding_method':
                setattr(params, key, value)
        
        # Encode
        encoder = Encode(params, logger=None)
        dna_codewords, info = encoder.encode(test_data)
        
        if not dna_codewords:
            raise ValueError(f"{encoding_method}: encode returned empty codewords")
        
        # Calculate cost
        total_cost = sum(len(seq) for seq in dna_codewords)
        
        # Cache result
        self._encode_cache[cache_key] = total_cost
        
        return total_cost
    
    @staticmethod
    def _make_cache_key(encoding_method: str, params_dict: Dict[str, Any]) -> str:
        """Create a cache key from encoding method and parameters."""
        # Sort params for consistent keys
        param_str = '|'.join(
            f"{k}={v}" for k, v in sorted(params_dict.items())
            if k != 'encoding_method'
        )
        return f"{encoding_method}:{param_str}"
    
    def clear_cache(self):
        """Clear encoding cache."""
        self._encode_cache.clear()


# Encoding-specific tuning strategies
class PerEncodingTuningStrategy:
    """Defines tuning strategies for each encoding."""
    
    @staticmethod
    def church() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for Church encoding."""
        return {
            'sequence_length': (100, 500),
            'rs_num': (0, 3),  # 0-3 redundant strands
        }
    
    @staticmethod
    def gcplus() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for GCPlus encoding."""
        return {
            'gcplus_k': (100, 250),  # Lower k = more oligos
        }
    
    @staticmethod
    def goldman() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for Goldman encoding."""
        return {
            'sequence_length': (100, 500),
        }
    
    @staticmethod
    def hedges() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for Hedges encoding."""
        return {
            'sequence_length': (100, 500),
            'hedges_coderate': (2, 8),  # Higher = more efficient
        }
    
    @staticmethod
    def max_density() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for MaxDensity encoding."""
        return {
            'codeword_length': (100, 1000),
        }
    
    @staticmethod
    def no_homopolymer() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for NoHomoPoly encoding."""
        return {
            'codeword_length': (100, 2000),
        }
    
    @staticmethod
    def wukong() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for Wukong encoding."""
        return {
            'sequence_length': (100, 500),
            'rs_num': (0, 3),
        }
    
    @staticmethod
    def yinyang() -> Dict[str, Tuple[float, float]]:
        """Parameters to tune for YinYang encoding."""
        return {
            'sequence_length': (50, 300),
            'max_content': (0.4, 0.7),  # GC content constraint
        }
    
    @classmethod
    def get_tuning_params(cls, encoding_method: str) -> Dict[str, Tuple[float, float]]:
        """Get tuning parameters for an encoding."""
        method = getattr(cls, encoding_method, None)
        if method is None:
            raise ValueError(f"No tuning strategy for {encoding_method}")
        return method()
