# Cost-Normalized Encoding Comparison Framework

## Overview

This framework allows you to compare DNA encodings based on **error resistance at equal cost**. The key insight is that different encodings produce different numbers of strands and strand lengths for the same input data. By normalizing the cost (total DNA length), you can fairly compare their ability to survive errors.

## Key Concepts

### Cost Definition

Cost is measured as **total DNA length** = sum of lengths of all DNA strands produced.

For example:
- Encoding A produces 100 strands of 300 bp each → cost = 30,000 bp
- Encoding B produces 200 strands of 150 bp each → cost = 30,000 bp
- Both encodings have the same cost, but different strand structures

### Workflow

1. **Cost Analysis Phase**
   - Encode a test dataset with each encoding using default parameters
   - Measure the total DNA length produced by each
   - Select a reference encoding (default: Goldman)

2. **Parameter Tuning Phase**
   - For each non-reference encoding, adjust its parameters to produce the same total DNA length as the reference
   - This is done by tuning the encoding's main cost-controlling parameter (e.g., `sequence_length`, `codeword_length`, etc.)

3. **Simulation Phase**
   - Run the full DNAbyte simulation pipeline for each encoding
   - Subject each to the same error conditions (synthesis, storage, sequencing)
   - Measure recovery success rate (can original data be decoded after errors?)

4. **Analysis**
   - Compare error resistance across encodings at equal cost
   - Identify which encodings are most robust given a fixed DNA synthesis budget

## File Structure

```
simulations/
├── cost_comparison.py              # Core framework classes
├── run_cost_comparison.py           # Example usage script
├── cost_normalized_results/         # Output directory
│   ├── cost_analysis.json           # Cost metrics for each encoding
│   ├── simulation_results.json      # Simulation results
│   └── comparison_YYYYMMDD_HHMMSS.log
└── ...
```

## Usage

### Basic Usage

```python
from simulations.run_cost_comparison import CostNormalizedComparison

# Run comparison
comparison = CostNormalizedComparison()
cost_results, sim_results = comparison.run(
    encodings=['church', 'gcplus', 'goldman', 'hedges', 'max_density', 
               'no_homopolymer', 'wukong', 'yinyang'],
    test_file='./tests/testfiles/test_data.bin',
    num_runs=3,
    reference_encoding='goldman'
)

# Results are saved to ./simulations/cost_normalized_results/
```

### Advanced: Customizing Parameters

```python
from simulations.cost_comparison import CostComparisonRunner

base_params = {
    'sequence_length': 250,           # Override default
    'max_homopolymer': 4,
    'synthesis_method': 'error_prone',
    'storage_conditions': {'temp': 4, 'time_days': 100},
    'sequencing_method': 'nanopore',
}

cost_runner = CostComparisonRunner()
results = cost_runner.run_cost_analysis(
    encodings=['church', 'goldman', 'wukong'],
    test_data=test_data,
    base_params=base_params,
    reference_encoding='goldman'
)
```

### Advanced: Implementing Custom Tuning

To implement smarter parameter tuning for an encoding, modify the corresponding `_tune_*` method in `ParameterTuner`:

```python
def _tune_hedges(self, target: int, params: Dict, test_data: str) -> Dict:
    """Tune Hedges to achieve target total DNA length."""
    
    # Binary search for the right sequence_length
    lo, hi = 100, 1000
    best_params = params.copy()
    best_cost = float('inf')
    
    while lo <= hi:
        mid = (lo + hi) // 2
        test_params = params.copy()
        test_params['sequence_length'] = mid
        
        # Would need to encode here to get actual cost
        # Simplified for now
        cost = self._estimate_cost(test_params, test_data)
        
        if cost < target:
            best_params = test_params
            best_cost = cost
            lo = mid + 1
        else:
            hi = mid - 1
    
    return best_params

def _estimate_cost(self, params: Dict, test_data: str) -> int:
    """Estimate the cost without actually encoding (faster)."""
    # Implement encoding-specific cost estimation
    pass
```

## Encoding-Specific Tuning Parameters

Here's a guide for implementing full parameter tuning for each encoding:

### 1. Church
- **Main tuner**: `sequence_length` (default 200)
- **Secondary tuners**: `rs_num` (redundancy), `add_redundancy`
- **Relationship**: Longer sequences → fewer strands → lower cost
- **Constraints**: Must be compatible with input data size

### 2. GCPlus
- **Main tuner**: `gcplus_k` (default 168) — information bits per oligo
- **Relationship**: Lower k → more oligos needed → higher cost
- **Constraints**: k affects oligo length due to error correction overhead

### 3. Goldman
- **Main tuner**: `sequence_length` (default 200)
- **Note**: Simplest encoding, fixed 4x overlapping overhead
- **Relationship**: Longer sequences → fewer strands → lower cost

### 4. Hedges
- **Main tuner**: `sequence_length` (default 300)
- **Secondary tuner**: `hedges_coderate` (default 3)
- **Relationship**: Longer sequences OR higher code rate → lower cost
- **Note**: Code rate is efficiency factor; higher = more efficient

### 5. MaxDensity
- **Main tuner**: `codeword_length` (default 500)
- **Relationship**: Longer codewords → fewer codewords → lower cost
- **Note**: Uses 2 bits/nucleotide, near theoretical maximum
- **Formula**: cost ≈ data_bits / 2 + barcode_overhead

### 6. NoHomoPoly
- **Main tuner**: `codeword_length` (default 500)
- **Secondary tuners**: `dna_barcode_length`, error correction parameters
- **Relationship**: Longer codewords → lower cost
- **Note**: Odd-even encoding adds overhead vs. MaxDensity

### 7. Wukong
- **Main tuner**: `sequence_length` (default 200)
- **Secondary tuners**: GC constraints (`min_gc`, `max_gc`)
- **Relationship**: Longer sequences OR relaxed GC constraints → lower cost
- **Note**: Tighter GC constraints reduce encoding efficiency

### 8. YinYang
- **Main tuner**: `sequence_length` (default 120)
- **Secondary tuner**: `max_content` (GC constraint, default 0.6)
- **Relationship**: Longer sequences OR relaxed GC constraint → lower cost
- **Note**: Search count is computational, doesn't affect DNA cost

## Expected Output

### cost_analysis.json
```json
{
  "church": {
    "num_strands": 156,
    "avg_strand_length": 223.5,
    "total_dna_length": 34866,
    "parameters": {...}
  },
  "goldman": {
    "num_strands": 200,
    "avg_strand_length": 174.3,
    "total_dna_length": 34860,
    "parameters": {...}
  },
  ...
}
```

### comparison log output
```
Encoding Costs:
  church              : strands=  156, avg_len= 223.5, total_dna=  34866
  gcplus              : strands=  234, avg_len= 149.0, total_dna=  34866
  goldman             : strands=  200, avg_len= 174.3, total_dna=  34860
  ...

Error Resistance (from simulation):
  church              :  85.3% success (17/20 runs)
  gcplus              :  92.1% success (18/20 runs)
  goldman             :  78.9% success (15/20 runs)
  ...
```

## Next Steps

1. **Implement full parameter tuning** for each encoding
   - Create `_tune_*` methods that use actual encoding calls for accurate cost estimation
   - Use binary search or gradient descent to find optimal parameters
   - Cache encoding results to avoid repeated calls

2. **Add error model variations**
   - Test with different error rates (low, medium, high)
   - Test with different error types (substitution, insertion, deletion)
   - Create trade-off curves: cost vs. error resistance

3. **Analyze results**
   - Which encoding is most robust at low cost?
   - Which is most efficient at high reliability?
   - Are there encodings that dominate others (better resistance at lower cost)?

4. **Optimize further**
   - Combine best parameters from multiple encodings
   - Implement adaptive parameter selection based on error model

## Implementation Tips

- **Performance**: Encoding can be slow. Cache results aggressively.
- **Accuracy**: Real cost depends on actual encoding, not estimates. Use actual encode() calls in parameter tuning.
- **Debugging**: Enable debug logs in simulation to see detailed step-by-step results.
- **Testing**: Start with a small test file (~1KB) to quickly iterate on parameters.

## References

- See [DNAbyte documentation](../docs/) for encoding details
- See [Simulation parameters](./simulation_parameters.py) for full list of tunable parameters
