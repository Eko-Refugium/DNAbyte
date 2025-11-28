# Pipeline Steps

## Step 0: Parameter Configuration (`Params` class)

The `Params` class encapsulates all configuration parameters required to simulate a full DNA data storage workflow. The parameters can be grouped into the steps of the pipeline they concern:

### Parameter Groups

- **File & Identification**: `name`, `seed`, `debug`
- **Encoding**: `encoding_scheme`, `assembly_structure`, `library_name`, `codeword_length`, `sigma_amount`, `percent_of_symbols`
- **Error Correction**: `inner_error_correction`, `outer_error_correction`, `reed_solo_percentage`, `ltcode_header`, `index_carry_length`
- **Synthesis**: `synthesis_method`, `mean`, `vol`, `std_dev`, `hybridisation_steps`
- **Storage**: `years`, `storage_conditions`
- **Sequencing**: `sequencing_method`, `iid_error_rate`

More details about the parameters will be provided in the respective sections. 

### Example Usage

```python
from dnabyte.params import Params

params = Params(
    name="test1",
    encoding_scheme="linear_encoding",
    assembly_structure="linear_assembly",
    outer_error_correction="reedsolomon",
    reed_solo_percentage=10,
    dna_barcode_length=6,
    codeword_length=100,
    library_name="standard_library.txt"
)
```

## Step 1: Binarize and Compress Data

To initiate the binarization process, instantiate a [`DNADS`](reference.md#dnads) object with a list of file paths. Then use this object to create a [`RawData`](reference.md#rawdata) object.

```python
from dnabyte.data import DNADS, RawData

data_obj = DNADS(["mi_dna_disc_logo.png"])
data_raw = RawData(data_obj)
```

The [`RawData`](reference.md#rawdata) class compresses all input files into a temporary `.tar.gz` archive, reads it as a binary stream, and converts it into a bitstring. This binary string is stored in `data_raw.data`, with its length available via `data_raw.length`. The temporary archive is deleted after use. The [`RawData`](reference.md#rawdata) class is a subclass of [`DNADS`](reference.md#dnads). It represents the compressed binary form of one or more input files. It is used to initiate the encoding process.

### Alternative Instantiation from Bitstream

For testing purposes, a [`RawData`](reference.md#rawdata) object can also be created from a binary string:

```python
raw_data = RawData.from_bitstream('100001010101')
```
This method is especially useful in simulations or unit tests where reading actual files is not required.

---

## Step 2: Encode Data

The process of encoding raw data involves two key steps:

* Instantiate an [`Encode`](reference.md#Encode') Object by providing all the required parameters, stored in the ['Params'](reference.md#params).

* Apply the ['encode'](reference.md#encode) Method:
Use the encode method on a RawData object to transform the raw binary data into DNA sequences. This method applies the specified encoding strategy and optionally incorporates error correction mechanisms.

This two-step approach is a consistent design pattern used throughout the framework when transitioning between different data classes.

### Usage

```python
from dnabyte.encode import Encode

enc = Encode(params=params, logger=self.simlogger)
data_enc, info = enc.encode(data_raw)
```

### Encoding Parameters

| Assembly Structure      | Encoding Scheme                   | Encoder Used               |
|-------------------------|------------------------------------|----------------------------|
| `linear_assembly`       | `linear_encoding`                  | `EncodeLinearChain`        |
| `linear_assembly`       | `binomial_encoding`                | `EncodeLinearBinom`        |
| `positional_assembly`   | `linear_encoding`                  | `EncodePolyChain`          |
| `positional_assembly`   | `binomial_encoding`                | `EncodePolyBinom`          |
| `synthesis`             | `max_density_encoding`             | `EncodeMaxDensity`         |
| `synthesis`             | `no_homopolymeroddeven_encoding`   | `EncodeNoHomoPolyOddEven`  |

### Error Correction Parameters

| Parameter              | Description                                                                 |
|------------------------|-----------------------------------------------------------------------------|
| `inner_error_correction` | Selects inner-layer error correction (e.g., `ltcode`, `None`)               |
| `ltcode_header`          | Number of bits used for Luby Transform code headers |
| `percent_of_symbols`     | Percentage of symbols included in each encoded LT packet                    |
| `index_carry_length`     | Number of bits used to carry the LT-code packet index                       |
| `outer_error_correction` | Selects outer-layer error correction (e.g., `reedsolomon`, `None`)          |
| `reed_solo_percentage`   | Redundancy percentage for Reed-Solomon coding                               |

These error correction methods are commonly used for synthesis-based methods, where the information lies in the sequence of nucleotides. They can, however, also be used for assembly based methods, where the information lies in the sequence of motives (short nucleotide sequences). Bare in mind, that assembly-based methods exhibit an implicit error correction that results from the fact that it lies in a restricted sequence space and that in most cases error correction is not required.  

## Returns

- An `EncodedData` object containing the encoded oligonucleotide sequences.
- An `info` dictionary with encoding metrics and metadata.

---

## Step 3: Simulate Synthesis

**Class Involved**: `SimulateAssembly`

Simulates the synthesis process of converting digital data into DNA oligonucleotides. Essentially, there are two different snythesis strategies that can be applied, which are determined by the parameter 'assembly_structure'. The option `synthesis` describes the nucleotide-wise chemical synthesis, which is based on the methods from Schwarz et al. The options `linear_assembly` or `positional_assembly` describe the synthesis by ligating short oligos from a so called library. This approach is based on a simulation of the hybridisation and ligation of the molecules. 

### Usage

```python
from dnabyte.assembly import SimulateAssembly

syn = SimulateAssembly(params, logger=self.simlogger)
data_syn, info = syn.simulate(data_enc)
```

### Synthesis Parameters

| Synthesis Approach              | Method                                        | Method ID  |
|---------------------------------|-----------------------------------------------|------------|
| Column synthesized oligos       |                                               |            |
|                                 | MutS                                          | `68`       |
|                                 | Consensus shuffle                             | `69`       |
|                                 | ErrASE                                        | `6`        |
|                                 | No error correction                           | `71`       |
| Microarray based oligo pools    |                                               |            |
|                                 | Oligo hybridization based error correction    | `4`        |
|                                 | High-temperature ligation/hybridization based | `5`        |
|                                 | ErrASE                                        | `3`        |
|                                 | Nuclease-based                                | `7`        |
|                                 | NGS-based                                     | `70`       |
|                                 | No error correction                           | `71`       |

For assembly based approaches only the following parameters must be given:

| Parameter        | Description                                                        | 
|------------------|--------------------------------------------------------------------|
| `mean`           | the mean number of copies of an oligo                              |
| `std_dev`        | the standard deviation of copies of an oligo                       | 
| `vol`            | the volume of the solution                                         | 
| `hybridisation`  | the number of random hybridisation steps (proxy for reaction time) |

### Reference

[1] Schwarz, Michael, et al. "MESA: automated assessment of synthetic DNA fragments and simulation of DNA synthesis, storage, sequencing and PCR errors." *Bioinformatics* 36.11 (2020): 3322-3326.


## Step 4: Simulate Storage

The `SimulateStorage` class models the long-term degradation of DNA sequences stored under various environmental conditions. DNA degradation is modeled probabilistically based on empirical decay rates under three conditions:

- **Room temperature**: approximated half-life of 521 years (decay rate ≈ 0.00133)
- **Permafrost**: 5.5×10⁻⁶ per nucleotide per year  
  *Source: Allentoft et al., 2012*
- **Biogene capsule**: 1×10⁻⁷ per nucleotide per year  
  *Source: Coudy et al., 2021*

### Usage

```python
from dnabyte.storage import SimulateStorage

sto = SimulateStorage(params, logger=self.simlogger)
data_sto, info = sto.simulate(data_syn)
```
### Storage Parameters

| Parameter             | Type     | Description                                                                 |
|-----------------------|----------|-----------------------------------------------------------------------------|
| `years`               | `int`    | Number of years to simulate storage. Use `0` to disable simulation.        |
| `storage_conditions`  | `str`    | One of `'permafrost'`, `'room_temperature'`, `'biogene'`, or `'random'`.   |

The synthetic `'random'` mode allows base substitutions to simulate artificial mutation rates.

### Return Values

- `StoredData` object (filtered or mutated sequences)
- `info` dictionary containing:
  - `number_of_strand_breaks`: the number of sequences simulated to be degraded


### References

[1] Allentoft, M. E. et al. *The half-life of DNA in bone: measuring decay kinetics in 158 dated fossils.* Proc. R. Soc. B 279, 4724–4733 (2012).  
[2] Coudy, D. et al. *Long term conservation of DNA at ambient temperature: Implications for DNA data storage.* PLoS One 16.11 (2021): e0259868.

---


## Step 5: Simulate Sequencing

The `SimulateSequencing` class models the sequencing process by applying sequencing-specific error profiles to stored DNA oligonucleotides. This step reflects how technologies like Illumina or Nanopore introduce errors into DNA reads.



### Usage

```python
from dnabyte.sequencing import SimulateSequencing

simulator = SimulateSequencing(params)
sequenced_data, info = simulator.simulate(stored_data)
```


## Supported Modes

| Sequencing Method | Description                                                                 |
|-------------------|-----------------------------------------------------------------------------|
| `None`            | No sequencing simulation is performed.                                      |
| `'iid'`           | Applies an independent and identically distributed error model (custom rate). |
| `35` | Illumina: Single End |
| `36` | Illumina: Paired End |
| `39` | Nanopore: 1D |
| `40` | Nanopore: 2D |
| `37` | PacBio: Subread|
| `38` | PacBio: CCS |


---
