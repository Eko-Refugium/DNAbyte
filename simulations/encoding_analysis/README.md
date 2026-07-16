# Encoding Analysis Module

Pre-simulation characterization of encoding properties to support parameter tuning.

## Purpose

Before running benchmarks, analyze each encoding to calculate metrics that DEPEND ON simulation parameters:
- **Data density**: bits per DNA base
- **Synthesis cost**: number of oligos required per unit data
- **Coverage**: average copy number of oligos
- **Error correction capability**: max tolerable error rate (depends on synthesis_method, mean, years, assembly_probability)
- **Assembly requirements**: reliability needed for synthesis success

## Key Insight

Error correction capability is NOT a fixed property of an encoding - it depends on:
- **Synthesis method**: MESA increases redundancy, assembly reduces it
- **Copy numbers**: Higher mean/std_dev improves tolerance
- **Storage duration**: Longer storage (years) reduces effectiveness
- **Assembly probability**: Lower probability reduces recovery success

## Workflow

1. **Run encoding-specific analyzer** → get metrics BASED ON your parameters
2. **Review breakdown** → understand which factors limit your capability
3. **Adjust parameters** → increase copies, use MESA synthesis, reduce storage time
4. **Re-analyze** → verify improvements
5. **Run benchmarks** with tuned parameters

## Usage Example

```python
from simulations.encoding_analysis import AnalyzeEncoding
from dnabyte.params import Params

# Analyze with specific parameters
params = Params(
    encoding_method='max_density',
    synthesis_method='mesa',        # Affects error capability!
    mean=20,                        # Affects error capability!
    years=1000,                     # Affects error capability!
    assembly_probability=0.90,      # Affects error capability!
)

analyzer = AnalyzeEncoding(params)
metrics = analyzer.calculate_all_metrics()
analyzer.print_report()  # Shows breakdown of all factors

recommendations = analyzer.recommend_parameters()
```

## Report Breakdown Example

```
Max tolerable error rate: 5.9%
Breakdown:
  Base resilience: 0.05           (wukong is less robust)
  × Synthesis adjustment: 1.30    (MESA provides redundancy)
  × Copy redundancy: 1.20         (20 mean copies helps)
  × Storage degradation: 0.80     (1000 years is challenging)
  × Assembly reliability: 0.95    (90% success rate is good)
  = Final estimate: 5.9%
```

Each factor is PARAMETER-DEPENDENT and transparent.

## Files

- `README.md` - This documentation
- `analyze_metrics.py` - Main AnalyzeEncoding class
- `example_usage.py` - Full workflow example showing parameter tuning

