# Simulation Framework

The `Simulation` class orchestrates the full end-to-end process of DNA-based data storage simulation in DNAbyte. It simulates encoding, synthesis, storage, sequencing, error correction, and decoding of digital files using oligonucleotide sequences.

## Overview

The simulation pipeline consists of the following steps:

1. Binarize input data into a compressed bitstream
2. Encode the bitstream into oligonucleotide codewords
3. Simulate synthesis errors
4. Simulate storage degradation
5. Simulate sequencing
6. Perform error correction
7. Decode data
8. Compare with original
9. Restore output files

## Class Definition

```python
class Simulation():
    def __init__(self, simulation_parameters, debug=False):
        ...
```

### Initialization

- A `Params` object defines the settings for each simulation.
- A logger records each step in a timestamped file.
- A unique job identifier (`self.job_identifier`) is assigned to each simulation run.

## Method: run()

```python
def run(self, paralel=False):
    ...
```

This method executes the entire simulation. Currently, parallel processing is not implemented but planned.

### Key Simulation Steps

#### Step 1: Binarize Data

```python
data_obj = DNADS([params.filename])
data_raw = RawData(data_obj)
```

Input files are compressed and transformed into a binary stream.

#### Step 2: Encode Data

```python
enc = Encode(params, logger=self.simlogger)
data_enc, info = enc.encode(data_raw)
```

Encoding schemes and barcoding are applied to generate codewords.

#### Step 3: Simulate Synthesis

```python
syn = SimulateAssembly(params, logger=self.simlogger)
data_syn, info = syn.simulate(data_enc)
```

Simulates oligo synthesis errors such as dropout and substitutions.

#### Step 4: Simulate Storage

```python
sto = SimulateStorage(params, logger=self.simlogger)
data_sto, info = sto.simulate(data_syn)
```

Models degradation over time such as strand breaks.

#### Step 5: Simulate Sequencing

```python
seq = SimulateSequencing(params, logger=self.simlogger)
data_seq, info = seq.simulate(data_sto)
data_seq = SequencedData(data_seq.data)
```

Sequencing errors are simulated on degraded oligos.

#### Step 6: Process Data

```python
data_cor, info = enc.process(data_seq)
```

Applies decoding-preprocessing and error correction.

#### Step 7: Decode Data

```python
data_dec, valid, info = enc.decode(data_cor)
```

Attempts to recover the original bitstream.

#### Step 8: Compare Results

```python
comparison, res = compare(data_dec, data_raw, logger=self.simlogger)
```

Checks if the decoded output matches the original data.

#### Step 9: Restore File

```python
RestoredData(data_dec, output_folder='./simulations/simfiles/', job_identifier=self.job_identifier)
```

Generates the output file from decoded data.

---

## Logging

Each simulation logs its steps to:

```
simulations/simlogs/job_<timestamp>.log
```

Results are also saved as:

```
simulations/simlogs/res_<timestamp>.pickle
```

---

## Output

A dictionary `results` is returned and contains:

- Step-wise performance info
- Final simulation status
- Duration, error counts, and decoding validation

---

## Planned Features

- Parallel processing support
- More complex storage/sequencing models
- CLI integration

---

## Usage Example

```python
from dnabyte.simulation import Simulation
sim = Simulation([params1, params2])
results = sim.run()
```

---

## Contact

For questions, issues, or contributions, open a GitHub issue or contact the maintainers.

