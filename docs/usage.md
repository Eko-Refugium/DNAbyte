# Pipeline Steps

## Step 0: Parameter Configuration (`Params` class)

The `Params` class encapsulates all configuration parameters required to simulate a full DNA data storage workflow. The parameters can be grouped into the steps of the pipeline they concern:

### Parameter Groups

- **File & Identification**: `name`, `file_paths`, `seed`, `debug`
- **Binarization**: `binarization_method`, `compression_level`
- **Encoding**: `encoding_method`, `assembly_structure`, `library_name`, `codeword_length`, `dna_barcode_length`, `codeword_maxlength_positions`, `percent_of_symbols`, `sigma_amount`
- **Error Correction**: `inner_error_correction`, `outer_error_correction`, `reed_solo_percentage`, `ltcode_header`, `index_carry_length`
- **Synthesis**: `synthesis_method`, `mesa_synthesis_id`, `mean`, `vol`, `std_dev`, `hybridisation_steps`
- **Storage**: `years`, `storage_conditions`
- **Sequencing**: `sequencing_method`, `sequencing_mesa_id`, `iid_error_rate`

More details about the parameters will be provided in the respective sections. 

### Example Usage

```python
from dnabyte.params import Params

params = Params(
    name="test1",
    file_paths=['my_file.txt'],
    binarization_method='compressed',
    encoding_method='max_density',
    outer_error_correction='reedsolomon',
    reed_solo_percentage=0.8,
    dna_barcode_length=34,
    codeword_length=501
)
```

### Creating Parameter Ranges

Create multiple parameter sets by varying one parameter:

```python
params_list = Params.params_range(
    name='storage_duration_test',
    file_paths=['test.txt'],
    encoding_method='max_density',
    binarization_method='compressed',
    storage_conditions='biogene',
    years=[1, 10, 100, 1000]  # Creates 4 Params objects
)
# Returns a list of 4 Params objects, one for each year value
```

---

## Step 1: Binarize and Compress Data

To initiate the binarization process, create a `Data` object with a list of file paths, then use `Binarize` to convert it to binary.

```python
from dnabyte.data_classes.base import Data
from dnabyte.binarize import Binarize

# Create data object
data_obj = Data(file_paths=["mi_dna_disc_logo.png"])

# Configure and run binarization
params = Params(binarization_method='compressed', compression_level=9)
binarizer = Binarize(params)
binary_code = binarizer.binarize(data_obj)

print(f"Binary length: {len(binary_code.data)} bits")
print(f"Binary data: {binary_code.data[:50]}...")
```

The `Binarize` class compresses all input files into a temporary `.tar.gz` archive (when using `compressed` method), reads it as a binary stream, and converts it into a bitstring. This binary string is stored in `binary_code.data`. The temporary archive is deleted after use.

### Binarization Methods

- **`'compressed'`**: Tar+gzip compression (recommended, requires `compression_level` parameter)
- **`'text'`**: Plain text encoding
- **`'default'`**: No compression

### BinaryCode Class

The result is a `BinaryCode` object containing:
- `data` (str): Binary string (e.g., "10010101...")
- `file_paths` (list): Original file paths

---

## Step 2: Encode Data

The encoding process transforms binary data into DNA sequences with optional error correction.

### Usage

```python
from dnabyte.encode import Encode

encoder = Encode(params)
nucleobase_code, info = encoder.encode(binary_code)

print(f"Number of codewords: {len(nucleobase_code.data)}")
print(f"Barcode length: {info['barcode_length']}")
```

### Encoding Methods

| Encoding Method | Assembly Structure | Description |
|----------------|-------------------|-------------|
| `max_density` | `synthesis` | Maximum information density encoding |
| `no_homopolymer` | `synthesis` | Homopolymer-free encoding |
| `linear_chain` | `linear_assembly` | Linear chain assembly |
| `linear_binom` | `linear_assembly` | Binomial linear assembly |
| `poly_chain` | `polymerase_assembly` | Polymerase chain assembly |
| `poly_binom` | `polymerase_assembly` | Binomial polymerase assembly |

### Error Correction Parameters

| Parameter              | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `inner_error_correction` | Selects inner-layer error correction: `'ltcode'` or `None`               |
| `ltcode_header`          | Number of bits used for Luby Transform code headers |
| `percent_of_symbols`     | Percentage of symbols included in each encoded LT packet                    |
| `index_carry_length`     | Number of bits used to carry the LT-code packet index                       |
| `outer_error_correction` | Selects outer-layer error correction: `'reedsolomon'` or `None`          |
| `reed_solo_percentage`   | Redundancy percentage for Reed-Solomon coding (0.0-1.0)                               |

Example with error correction:

```python
params = Params(
    file_paths=['data.txt'],
    encoding_method='max_density',
    binarization_method='compressed',
    codeword_length=501,
    dna_barcode_length=34,
    inner_error_correction='ltcode',
    ltcode_header=34,
    percent_of_symbols=2,
    index_carry_length=34,
    outer_error_correction='reedsolomon',
    reed_solo_percentage=0.8  # 80% redundancy
)
```

These error correction methods are commonly used for synthesis-based methods, where the information lies in the sequence of nucleotides. They can also be used for assembly-based methods, where the information lies in the sequence of motifs (short nucleotide sequences). Note that assembly-based methods exhibit implicit error correction from their restricted sequence space.

### NucleobaseCode Class

The result is a `NucleobaseCode` object. The structure of `data` depends on the encoding method:

- **Synthesis methods** (`max_density`, `no_homopolymer`): `data` is a list of DNA strings
- **Assembly methods** (`linear_chain`, `poly_chain`, etc.): `data` is nested lists of library elements

---

## Step 3: Simulate Synthesis

Simulates the synthesis process of converting digital data into DNA oligonucleotides. There are two different synthesis strategies determined by the `synthesis_method` parameter:

- **`'mesa'`**: Nucleotide-wise chemical synthesis based on the MESA model
- **`'assembly'`**: Synthesis by ligating short oligos from a library (hybridization simulation) 

### Usage

```python
from dnabyte.synthesize import SimulateSynthesis

synthesizer = SimulateSynthesis(params)
synthesized_dna, info = synthesizer.simulate(nucleobase_code)

print(f"Number of sequences: {len(synthesized_dna.data)}")
```

### MESA Synthesis Parameters

| Synthesis Approach              | Method                                        | Method ID  |
|---------------------------------|-----------------------------------------------|------------|
| **Column synthesized oligos**   |                                               |            |
|                                 | MutS                                          | `68`       |
|                                 | Consensus shuffle                             | `69`       |
|                                 | ErrASE                                        | `6`        |
|                                 | No error correction                           | `71`       |
| **Microarray based oligo pools**|                                               |            |
|                                 | Oligo hybridization based error correction    | `4`        |
|                                 | High-temperature ligation/hybridization based | `5`        |
|                                 | ErrASE                                        | `3`        |
|                                 | Nuclease-based                                | `7`        |
|                                 | NGS-based                                     | `70`       |
|                                 | No error correction                           | `71`       |

### Assembly Synthesis Parameters

For assembly-based approaches (`synthesis_method='assembly'`), the following parameters must be provided:

| Parameter            | Description                                                        | 
|----------------------|--------------------------------------------------------------------|
| `mean`               | Mean number of copies of an oligonucleotide                        |
| `std_dev`            | Standard deviation of copy numbers                                 | 
| `vol`                | Volume of the solution (in liters)                                 | 
| `hybridisation_steps`| Number of random hybridization steps (proxy for reaction time)    |

Example:
```python
from scipy.constants import Avogadro

params = Params(
    synthesis_method='mesa',
    mesa_synthesis_id=68,
    mean=20,
    vol=1000000 / Avogadro,
    std_dev=1,
    hybridisation_steps=10000
)
```

### InSilicoDNA Class

The result is an `InSilicoDNA` object containing:
- `data` (list): List of DNA sequence strings
- Methods for statistics and string representation

### Reference

[1] Schwarz, Michael, et al. "MESA: automated assessment of synthetic DNA fragments and simulation of DNA synthesis, storage, sequencing and PCR errors." *Bioinformatics* 36.11 (2020): 3322-3326.

---

## Step 4: Simulate Storage

The `SimulateStorage` class models the long-term degradation of DNA sequences stored under various environmental conditions. DNA degradation is modeled probabilistically based on empirical decay rates under three conditions:

- **`'biogene'`**: 1×10⁻⁷ per nucleotide per year (Biogene capsule storage)
- **`'permafrost'`**: 5.5×10⁻⁶ per nucleotide per year  
- **`'roomtemperature'`**: approximated half-life of 521 years (decay rate ≈ 0.00133)
- **`'random'`**: Custom decay for testing

### Usage

```python
from dnabyte.store import SimulateStorage

storage = SimulateStorage(params)
stored_dna, info = storage.simulate(synthesized_dna)

print(f"Strand breaks: {info['number_of_strand_breaks']}")
```

### Storage Parameters

| Parameter             | Type     | Description                                                                 |
|-----------------------|----------|-----------------------------------------------------------------------------|
| `years`               | `int`    | Number of years to simulate storage. Use `0` to disable simulation.        |
| `storage_conditions`  | `str`    | One of `'permafrost'`, `'roomtemperature'`, `'biogene'`, or `'random'`.   |

The `'random'` mode allows base substitutions to simulate artificial mutation rates.

### Return Values

- `InSilicoDNA` object (filtered or mutated sequences)
- `info` dictionary containing:
  - `number_of_strand_breaks`: Number of sequences degraded


### References

[1] Allentoft, M. E. et al. *The half-life of DNA in bone: measuring decay kinetics in 158 dated fossils.* Proc. R. Soc. B 279, 4724–4733 (2012).  
[2] Coudy, D. et al. *Long term conservation of DNA at ambient temperature: Implications for DNA data storage.* PLoS One 16.11 (2021): e0259868.

---

## Step 5: Simulate Sequencing

The `SimulateSequencing` class models the sequencing process by applying sequencing-specific error profiles to stored DNA oligonucleotides.

### Usage

```python
from dnabyte.sequence import SimulateSequencing

sequencer = SimulateSequencing(params)
sequenced_dna, info = sequencer.simulate(stored_dna)

print(f"Errors introduced: {len(info)}")
```

### Sequencing Methods

| Sequencing Method | ID/Parameter | Description                                                                 |
|-------------------|--------------|-----------------------------------------------------------------------------|
| `None`            | -            | No sequencing simulation is performed.                                      |
| `'iid'`           | `iid_error_rate` | Independent and identically distributed error model (specify custom rate). |
| `'mesa'`          | `35` | Illumina: Single End |
| `'mesa'`          | `36` | Illumina: Paired End |
| `'mesa'`          | `39` | Nanopore: 1D |
| `'mesa'`          | `40` | Nanopore: 2D |
| `'mesa'`          | `37` | PacBio: Subread|
| `'mesa'`          | `38` | PacBio: CCS |

### Examples

```python
# MESA sequencing with Illumina
params = Params(
    sequencing_method='mesa',
    sequencing_mesa_id=35  # Illumina Single End
)

# IID sequencing with custom error rate
params = Params(
    sequencing_method='iid',
    iid_error_rate=0.01  # 1% error rate
)
```

### Return Values

- `InSilicoDNA` object with sequencing errors
- `info` list containing error information

---

## Step 6: Process Data

Align and prepare sequences for decoding. This step processes sequenced reads to prepare them for error correction and decoding.

### Usage

```python
# Process sequenced data
processed_code, info = encoder.process(sequenced_dna)

print(f"Processed {len(processed_code.data)} codewords")
```

**Note:** For assembly methods requiring sequence alignment (`linear_binom`, `poly_binom`), this step requires **MMseqs2** to be installed on your system.

### Return Values

- `NucleobaseCode` object ready for decoding
- `info` dictionary with processing metrics

---

## Step 7: Decode Data

Decode DNA sequences back to binary with error correction.

### Usage

```python
decoded_binary, valid, info = coder.decode(processed_code)

if valid:
    print("Decoding successful!")
    print(f"Recovered {len(decoded_binary.data)} bits")
else:
    print("Decoding failed - error correction could not recover data")
```

### Return Values

- `BinaryCode` object with recovered binary data
- `valid` (bool): `True` if error correction successfully recovered original data
- `info` dictionary with decoding metrics

The `valid` flag indicates whether the decoded data matches the original after error correction. If too many errors were introduced during synthesis/storage/sequencing, `valid` will be `False`.

---

## Step 8: Restore Files

Convert binary back to original files.

### Usage

```python
if valid:
    restored_data, info = binarizer.debinarize(decoded_binary)
    print(f"Files restored: {restored_data.file_paths}")
```

### Return Values

- `Data` object with restored files
- `info` dictionary with restoration metrics

The files are restored to their original format and can be accessed from the file system.

---

## Minimal Example (No Error Channels)

```python
from dnabyte.params import Params
from dnabyte.data_classes.base import Data
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode

params = Params(
    file_paths=['test.txt'],
    encoding_method='max_density',
    binarization_method='compressed',
    synthesis_method=None,
    storage_conditions=None,
    sequencing_method=None
)

data = Data(file_paths=['test.txt'])
binarizer = Binarize(params)
encoder = Encode(params)

binary = binarizer.binarize(data)
encoded, _ = encoder.encode(binary)
decoded, valid, _ = encoder.decode(encoded)

if valid:
    restored, _ = binarizer.debinarize(decoded)
    print("Success!")
```

---

## Additional References

[1] Schwarz, M., et al. "MESA: automated assessment of synthetic DNA fragments and simulation of DNA synthesis, storage, sequencing and PCR errors." *Bioinformatics* 36.11 (2020): 3322-3326.

[2] Allentoft, M.E., et al. "The half-life of DNA in bone: measuring decay kinetics in 158 dated fossils." *Proc. R. Soc. B* 279 (2012): 4724–4733.

[3] Coudy, D., et al. "Long term conservation of DNA at ambient temperature: Implications for DNA data storage." *PLoS One* 16.11 (2021): e0259868.
