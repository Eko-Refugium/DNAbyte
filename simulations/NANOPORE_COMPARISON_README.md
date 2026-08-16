# Nanopore Encoding Comparison Simulation

This simulation compares different DNA encoding methods' resilience to nanopore sequencing errors.

## Overview

The simulation:
1. Creates a random text file for encoding
2. Encodes it using different encoding methods (Goldman, Church, GC+, Max Density, No Homopolymer, Wukong, YinYang)
3. Simulates nanopore sequencing (both 1D and 2D methods)
4. Tracks insertions, deletions, and substitutions introduced during sequencing
5. Generates comparison bar charts

## Usage

### Basic Run

```bash
cd DNAbyte
python simulations/nanopore_encoding_comparison.py
```

### Configuration

Edit the script's `__main__` section to customize:

- **Encodings to test**: Comment out encodings you don't want to test
- **Nanopore methods**: `[39, 40]` where 39 = 1D, 40 = 2D
- **Number of runs**: Start with 3-5 for quick testing, increase to 10-20 for final results

Example:
```python
encodings = [
    'goldman',
    'church',
    # 'gcplus',      # Comment out to skip
    'max_density',
    # 'no_homopolymer',
    'wukong',
    'yinyang'
]
num_runs = 5  # Adjust as needed
```

## Output

Results are saved in `./simulations/nanopore_comparison/`:

### Files Generated

1. **comparison_results.json** - Raw results data
2. **nanopore_success_rates.png** - Success rate comparison
3. **nanopore_error_types.png** - Error types (insertions, deletions, substitutions)
4. **nanopore_detailed_comparison.png** - Side-by-side 1D vs 2D comparison

### Console Output

- Progress for each encoding and method
- Per-run error statistics (insertions, deletions, substitutions)
- Summary table with aggregated results
- Success rates for each encoding+method combination

## Understanding the Results

### Success Rate
Percentage of runs where the original data was successfully recovered after nanopore sequencing.

### Error Rates
Measured per 1000 bases:
- **Insertions**: Extra bases added during sequencing
- **Deletions**: Bases removed during sequencing
- **Substitutions**: Bases changed to different bases

### Nanopore Methods
- **1D**: Sequences one strand (higher error rate, ~20%)
- **2D**: Sequences both strands (lower error rate, ~13%)

## Expected Runtime

With default settings (7 encodings × 2 methods × 5 runs = 70 simulations):
- **Quick test** (500 byte file, 5 runs): ~10-15 minutes
- **Full test** (500 byte file, 20 runs): ~30-40 minutes

## Troubleshooting

### Missing Dependencies
```bash
pip install numpy matplotlib scipy tqdm
```

### File Not Found
Make sure you run from the `DNAbyte` directory:
```bash
cd DNAbyte
python simulations/nanopore_encoding_comparison.py
```

### Memory Issues
Reduce the number of runs or test fewer encodings at once.

## Customization

### Different File Size
Edit in the script:
```python
create_random_file(test_file_path, size_bytes=1000)  # Change from 500
```

### Different Sequencing Parameters
Modify `base_params` in the `run_comparison` method:
```python
base_params = {
    'sequence_length': 150,  # Change sequence length
    'mean': 15,             # Change coverage depth
    # ... other parameters
}
```

## Analysis Tips

1. **Most Resilient Encoding**: Highest success rate + lowest error counts
2. **Best for Nanopore**: Compare 1D vs 2D performance
3. **Error Type Analysis**: Some encodings may be better at preventing specific error types
4. **Trade-offs**: Consider both success rate AND total errors introduced

## References

- Mesa sequencing simulation: Based on MESA error model
- Nanopore error profiles: From Weirather et al., 2017
- Encoding methods: Various DNA storage encoding schemes
