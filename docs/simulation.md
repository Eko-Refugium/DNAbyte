# Simulation Framework

The `Simulation` class in `simulations/simulation.py` orchestrates automated end-to-end DNA storage simulations. It executes the complete pipeline from file binarization through synthesis, storage, sequencing, and decoding, while collecting detailed metrics and logging results.

---

## Overview

The simulation framework enables:

- **Batch Processing**: Run multiple parameter configurations in sequence
- **Automated Logging**: Comprehensive logs saved to `simulations/simlogs/`
- **Progress Tracking**: Real-time progress bars via tqdm
- **Error Handling**: Graceful failure handling with detailed error logs
- **Result Serialization**: Results saved as pickle files for later analysis
- **Reproducibility**: Optional random seed control for deterministic simulations

---

## Architecture

### Pipeline Steps

The simulation executes all steps of the DNAbyte data storage pipeline controlled by a params object:

Each step:
- Logs start time, status, duration, and key metrics
- Returns an `info` dictionary with step-specific results
- Handles errors gracefully without crashing the entire simulation
- Saves results for later analysis

---

## Basic Usage

### Simple Example

```python
from dnabyte.params import Params
from simulation import Simulation

# Configure simulation parameters
params = Params(
    name='my_simulation',
    file_paths=['data.txt'],
    encoding_method='max_density',
    binarization_method='compressed',
    synthesis_method='mesa',
    mesa_synthesis_id=68,
    storage_conditions='biogene',
    years=10,
    sequencing_method='mesa',
    sequencing_mesa_id=68
)

# Run simulation
sim = Simulation([params])
results = sim.run()

# Access results
print(results)
```

### Parameter Ranges

Run simulations across multiple parameter values:

```python
from dnabyte.params import Params
from simulation import Simulation

# Create parameter range (varies storage duration)
params_list = Params.params_range(
    name='storage_duration_test',
    file_paths=['test.txt'],
    encoding_method='max_density',
    binarization_method='compressed',
    storage_conditions='biogene',
    years=[1, 10, 100, 1000],  # Creates 4 parameter sets
    synthesis_method='mesa',
    mesa_synthesis_id=68,
    sequencing_method='iid',
    iid_error_rate=0.01
)

# Run all simulations
sim = Simulation(params_list)
results = sim.run()

# Analyze results
for name, data in results.items():
    if data['status'] == 'SUCCESS':
        print(f"{name}: Encoding took {data['step2']['duration']:.2f}s")
```

---

## Initialization

### Constructor

```python
Simulation(simulation_parameters, debug=False)
```

**Parameters:**
- `simulation_parameters` (list): List of `Params` objects or single `Params` object
- `debug` (bool): Enable debug mode (currently unused)

**Initialization actions:**
1. Creates logger instance
2. Generates unique job identifier (timestamp: `YYYYMMDD_HHMMSS`)
3. Creates log file: `simulations/simlogs/job_{timestamp}.log`
4. Sets up log formatting

---

## Running Simulations

### run() Method

```python
results = sim.run(paralel=False)
```

**Parameters:**
- `paralel` (bool): Enable parallel processing (not yet implemented)

**Returns:**
- `dict`: Nested dictionary with results for each parameter set

**Result Structure:**
```python
{
    'simulation_name_1': {
        'status': 'SUCCESS' or 'FAILURE',
        'step1': {'duration': 0.05, 'length_of_bitsream': 2144},
        'step2': {'duration': 0.67, 'barcode_length': 34, ...},
        'step3': {'duration': 0.12, ...},
        # ... more steps ...
    },
    'simulation_name_2': { ... }
}
```

---

## Logging

### Log File Structure

```
simulations/simlogs/job_20251211_163000.log
```

**Log Format:**
```
2025-12-11 16:30:00,123 - INFO - SIMULATION SETTING: my_simulation
2025-12-11 16:30:00,124 - INFO - Params:
  name: my_simulation
  file_paths: ['test.txt']
  ...
2025-12-11 16:30:00,125 - INFO - STEP01: BINARIZE DATA
2025-12-11 16:30:00,130 - INFO - STATUS: SUCCESS
2025-12-11 16:30:00,130 - INFO - DURATION: 0.00 seconds
```

### Log Levels

- **INFO**: Step execution, metrics, success messages
- **ERROR**: Failures, exceptions, error traces

---

## Result Serialization

Results are automatically saved as pickle files:

```
simulations/simlogs/res_20251211_163000.pickle
```

### Loading Results

```python
import pickle

with open('simulations/simlogs/res_20251211_163000.pickle', 'rb') as f:
    results = pickle.load(f)

# Analyze results
for name, data in results.items():
    if data['status'] == 'SUCCESS':
        total_time = sum(
            step['duration'] 
            for step in data.values() 
            if isinstance(step, dict) and 'duration' in step
        )
        print(f"{name}: Total time = {total_time:.2f}s")
```

---

## Error Handling

### Graceful Failures

Each step has independent error handling:

```python
try:
    # Execute step
    data_enc, info = coder.encode(binary_code)
except Exception as e:
    self.simlogger.info('STATUS: ERROR')
    self.simlogger.error('TYPE: %s', str(e))
    self.simlogger.error(traceback.format_exc())
    results[name]['status'] = 'FAILURE'
    continue  # Skip to next parameter set
```

**Behavior:**
- Individual step failures don't crash the simulation
- Error logged with full traceback
- Simulation continues with next parameter set
- Failed simulations marked with `status: 'FAILURE'`

---

## Progress Tracking

Real-time progress bar using tqdm:

```
Simulation Progress: 45%|████████      | 45/100 [02:15<02:45,  3.05param/s]
```

Shows:
- Percentage complete
- Current/total parameter sets
- Elapsed time
- Estimated time remaining
- Processing rate (params/second)

---

## Reproducibility

### Random Seed Control

```python
params = Params(
    name='reproducible_test',
    seed=42,  # Set random seed
    file_paths=['test.txt'],
    # ... other params ...
)

sim = Simulation([params])
results = sim.run()
```

**Behavior:**
- If `seed` is provided, sets `random.seed(params.seed + counter)`
- Ensures deterministic behavior for stochastic steps
- Different counter for each parameter set in batch

---

## Complete Example: Success Rate Analysis

```python
from dnabyte.params import Params
from simulation import Simulation

# Test error correction limits
params_list = Params.params_range(
    name='success_rate_analysis',
    file_paths=['test.txt'],
    encoding_method='max_density',
    binarization_method='compressed',
    outer_error_correction='reedsolomon',
    reed_solo_percentage=0.8,
    synthesis_method='mesa',
    mesa_synthesis_id=68,
    storage_conditions='biogene',
    years=[1, 10, 100, 1000, 10000],  # Vary storage duration
    sequencing_method='iid',
    iid_error_rate=0.01,
    seed=42
)

# Run 100 repeats for each year value
all_params = []
for _ in range(100):
    all_params.extend(params_list)

sim = Simulation(all_params)
results = sim.run()

# Calculate success rates
from collections import defaultdict
success_counts = defaultdict(lambda: {'success': 0, 'failure': 0})

for name, data in results.items():
    # Extract year from params (would need to track this)
    year = ...  # Extract from name or params
    if data['status'] == 'SUCCESS':
        success_counts[year]['success'] += 1
    else:
        success_counts[year]['failure'] += 1

# Plot success rates
for year, counts in success_counts.items():
    rate = counts['success'] / (counts['success'] + counts['failure'])
    print(f"Year {year}: Success rate = {rate:.2%}")
```

---

## Planned Features

- **Parallel Processing**: Multi-threaded/multi-process execution
- **Checkpointing**: Resume interrupted simulations
- **CLI Integration**: Command-line interface for simulations
- **Real-time Monitoring**: Web dashboard for live progress
- **Custom Callbacks**: User-defined functions at each step

---

## Tips and Best Practices

1. **Start Small**: Test with small files before large-scale simulations
2. **Use Seeds**: Set `seed` parameter for reproducible results
3. **Monitor Logs**: Check log files for detailed error information
4. **Batch Wisely**: Balance batch size with memory constraints
5. **Save Results**: Always check that pickle files are saved successfully
6. **Clean Up**: Remove test files after simulations complete

---

## Troubleshooting

### Common Issues

**"FileNotFoundError: No such file or directory"**
- Ensure file paths are relative to project root
- Check that input files exist before running

**"MMseqs2 not installed"**
- Required for `linear_binom` and `poly_binom` encodings
- Install: `brew install mmseqs2` (macOS) or equivalent

**"Decoding failed - error correction could not recover data"**
- Expected behavior when errors exceed correction capacity
- Increase `reed_solo_percentage` or reduce error rates

**Memory errors with large files**
- Process files in smaller batches
- Reduce `file_paths` list size

---

## Contact

For questions, issues, or contributions, open a GitHub issue or contact the maintainers.
