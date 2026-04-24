
## Data Classes

Data Classes are the stored intermediate states of the pipeline, this means between every step of the pipeline the structure of the ingoing and resulting data is checked to be compatible with the rest of the pipeline. Each data class represents a structure that arises in the pipeline. The type system prevents mismatches, ensuring you can't accidentally pass binary data to a function expecting DNA sequences.

Classes carry additional a dictionary that can be used by the modules to store metadata used in debugging and benchmarking.

### The Data throughout the pipeline

```
Data                   → Raw file content, filesystem representation
    ↓
BinaryCode            → Binary representation
    ↓
NucleobaseCode        → Encoded DNA sequences with structure needed in the synthesis/assembly
    ↓
InSilicoDNA           → Simulated physical DNA molecules
    ↓
InSilicoDNA           → After errors channels and sequencing
    ↓
NucleobaseCode        → Recovered encoded sequences
    ↓
BinaryCode            → Recovered binary representation
    ↓
Data                   → Restored files
```

### Core Data Classes Explained

#### **Data**
*Base class - raw digital files in native formats*

Pipeline entry/exit point representing filesystem content. Encapsulates file paths and calculates total size. All other data classes inherit from this base class, ensuring compatibility across the pipeline. Validates file existence and enforces size limits (default 1MB max per file).

#### **BinaryCode (extends Data)**
*Binary representation as bitstream*

Stores data as string of '0' and '1' characters. Inherits from Data to maintain file_paths and size information. Provides validation ensuring only binary characters and comparison operations. Created by binarization methods or directly from binary strings.

#### **NucleobaseCode (extends Data)**
*Encoded DNA sequences*

- De novo synthesized encodings
Assumes a list structure containing each codeword once and maintains metadata of the encoding.

- Enzymatic assembly  
Assumes as a Multi-level list structure representing encoded data. Top level contains codewords, lower levels contain motifs/nucleotides organized by pooling hierarchy. Supports random assembly, iterative assembly, comparison, slicing, and reverse complement operations. Maintains encoding scheme metadata and biological constraints.

#### **InSilicoDNA (extends Data)**
*Simulated DNA molecules as oligonucleotide strings*

Contains list of DNA sequences (strings of 'A', 'C', 'G', 'T'). Used throughout synthesis simulation (in the case of de novo synthesis), storage simulation, and sequencing. Same class represents both designed molecules and error-containing molecules after physical channel simulation.

### Modularity

The Type system enforces that plugins accept and return appropriate classes. An error channel simulator only accepts `InSilicoDNA` and returns `InSilicoDNA`, preventing incompatible combinations. Within interface constraints, plugins have complete freedom. One encoding scheme might produce long heavily error resistant codes and store extensive metadata in `EncodedData` while another is minimalist. Both work as long as they satisfy the interface contract.

Clear input-output relationships mean stages compose naturally. Any synthesis simulator chains with any storage simulator because output type matches input type, enabling arbitrary combinations without custom glue code. Further data class boundaries provide natural unit testing seams. Test an encoder by feeding mock `BinaryCode` and validating the resulting `EncodedData`, without running the full pipeline.

### Practical Example: Following Data Through the Pipeline

To make this concrete, imagine storing an image file:

1. **Data**: 10,000-byte JPEG file with path "/path/to/photo.jpg"
2. **BinaryCode**: 80,000-bit string ('0' and '1' chars) plus original file_paths from Data
3. **NucleobaseCode**: Multi-level structure with ~500 DNA codewords organized hierarchically, each 200 nt
4. **InSilicoDNA**: Flat list of ~500 DNA sequence strings, ~2% now contain synthesis errors and redundancy
5. **InSilicoDNA** (after sequencing): Same structure but sequences now include sequencing errors and coverage variations (some sequences read multiple times, others missing)
6. **NucleobaseCode**: Reconstructed hierarchical structure during decoding
7. **BinaryCode**: bit string recovered (ideally identical to step 2)
8. **Data**: Restored JPEG file

Each transformation maintains the Data inheritance, ensuring all classes have compatible base methods.

---

## Pipeline Architecture

The pipeline consists of individual **steps**, each representing a stage in the DNA data storage lifecycle:

```
Digital File → Binarize → Encode → Synthesis → Storage → Misc Errors → 
Sequencing → Process → Decode → Compare → Restore File
```

---

## Step-by-Step Breakdown

### **STEP 1: Binarize Data** 
*Transform digital files into binary representation*

**Conceptual Purpose**: Convert arbitrary file formats to standardized binary representation for uniform downstream processing.

Binarization produces the input for the `BinaryCode` data class containing pure binary representation plus reconstruction metadata, abstracting away file format complexity.

**Modularity**: Plugins are to be added to the `binarization/` folder. The read file is to be named `binarize.py` and has to implement `binarize()` and `debinarize()` methods. Additional auxiliary files may be added in the folder. Strategies may include compression, encryption, or format-specific optimizations. Rest of pipeline operates on standardized binary, independent of source formats.

---

### **STEP 2: Encode Data**
*Convert binary data into DNA sequences with error correction*

**Conceptual Purpose**: Map binary to quaternary DNA alphabet while adding error resilience or any other required additions.

This step provides the input for the `NucleobaseCode` data class.

**Modularity**: Plugins in added in the `encoding/` folder and implement `encode()`, `decode()` and `processing()` methods in a file that must be named `encode.py`. The reason why the encoding plugin contains all three is because these are intertwined. `encode()` turns the binary string given into the quaternary code and if needed the assembly structure. Additional error correction methods that are more general that can be used in encodings are stored in `error_correction`. These can be e.g. fountain codes.

---

### **STEP 3: Simulate Synthesis**
*Model errors introduced during DNA strand creation*

**Conceptual Purpose**: Model transition from theoretical sequences to imperfect physical molecules based on synthesis chemistry error profiles.

Synthesis produces `InSilicoDNA` containing flat list of DNA sequence strings with potential synthesis errors.

**Modularity**: Plugins in `synthesis/` folder range from uniform error rates to platform-specific position-dependent models. Different error types (deletions, substitutions, insertions) have different impacts on encoding robustness. This step can be skipped for pure algorithmic testing. The read in file must be named `synthesize.py`

---

### **STEP 4: Simulate Storage**
*Model DNA degradation during physical storage*

**Conceptual Purpose**: Model how molecular damage accumulates over time based on storage conditions.

Storage simulation operates on `InSilicoDNA`, adding strand breaks, base modifications, and molecular loss while tracking fragment length distributions and damage locations.

**Modularity**: Plugins in `storage/` folder simulate different degradation mechanisms and environmental conditions (room temperature aging, archival storage at -80°C, etc.). Reveals trade-offs between sequence length and fragmentation resistance. Can be skipped for near-term data access scenarios. The read in file must be named `store.py`

---

### **STEP 5: Simulate Miscellaneous Errors**
*Model additional error sources throughout the workflow*

**Conceptual Purpose**: Model handling and processing steps (amplification, purification, library preparation) that introduce secondary errors, bias, or data loss.

Operates on `InSilicoDNA`, adding amplification errors, contamination, and sample loss while maintaining the same data structure.

**Modularity**: Plugins in `misc_errors/` folder provide error models that don't fit other categories. Multiple plugins can be chained together to combine PCR errors, contamination, and amplification bias. Can be omitted for simplified analyses focusing on primary error sources. The read in file must be named `err.py`

---

### **STEP 6: Simulate Sequencing**
*Model errors from DNA reading/sequencing process*

**Conceptual Purpose**: Simulate the retrieval operation and the challenge that stored DNA is never observed perfectly.

Sequencing operates on `InSilicoDNA`, producing observed sequences with sequencing errors.

**Modularity**: Plugins in `sequencing/` folder model different technologies with characteristic error profiles, read length distributions (Illumina, Nanopore, PacBio). The read in file must be named `sequence.py`

---

### **STEP 7: Process Data**
*Prepare sequenced reads for decoding*

**Conceptual Purpose**: Bridge the gap between raw sequencing output and structured decoder input by transforming noisy reads into clean, organized data structures.

Processing operates on `InSilicoDNA`, applying quality filtering, clustering, and consensus building to refine sequence quality and reduce redundancy.

**Modularity**: Plugins in are read in with the encoding and should be the step for cleaning up the reads, from simple quality filtering to machine learning, graph-based assembly, or statistical models could be implemented here.

---

### **STEP 8: Decode Data**
*Convert DNA sequences back to binary using error correction*

**Conceptual Purpose**: Reconstruct original binary data from processed sequences by exploiting structured redundancy added during encoding.

Decoding reconstructs `NucleobaseCode` from `InSilicoDNA` sequences, then applies error correction to produce `BinaryCode` with validation flag.

**Modularity**: Each encoding plugin in `encoding/` folder provides its corresponding decoder bundled in the same module. Different approaches use different algorithms based on the encoding.

---

### **STEP 9: Restore Data**
*Reconstruct original files from binary data*

**Conceptual Purpose**: Transform binary data back into concrete, usable files that can be opened and used.

Restoration reverses Step 1, using metadata in `BinaryCode` to reconstruct `Data` with all files, directory structure, and attributes via `debinarize()` operation.

**Modularity**: The plugin that created `BinaryCode` in Step 1 inverts the transformation. If compressed, decompresses; if encrypted, decrypts; if optimized, reverses optimizations. Provides human-level validation beyond binary correctness.

---

## Using the Pipeline

Each step can be tested independently and Swapped with alternative implementations to try different encoding schemes or error models by changing plugin files

The pipeline now needs to be configured. This can be done through the `Params` object. this gives a dictionary of all the parameters needed for the modules. For this reason every module needs a function called `attributes()` this takes the Params object and test if all necessary parameters are given for all modules that is loaded in. If not they can add default values this function in up to the user.

Examples of what the `Params` object can configure:
- **File inputs**: Specify which files to store, supporting single or multiple file scenarios
- **Encoding schemes**: Select from available encoding methods (Reed-Solomon, fountain codes, etc.)
- **Error correction methods**: Choose redundancy levels and error correction algorithms
- **Error rates and types**: Set synthesis/sequencing error probabilities and distributions
- **Storage conditions**: Define temperature, humidity, duration affecting degradation rates
- **Sequencing parameters**: Specify platform type, coverage depth, and quality thresholds

This enables systematic exploration of the design space to identify optimal configurations for specific applications.

## Modular Plugin Architecture

### Extensibility Through Folder Structure

Plugin-based architecture enables extending functionality through designated folders without modifying core code.

### How It Works

Pipeline automatically discovers and loads implementations from designated folders at runtime:

```
dnabyte/
├── binarization/          # Custom binarization methods
├── encoding/              # Encoding schemes
├── error_correction/      # Error correction algorithms
├── synthesis/             # Synthesis simulation models
├── storage/               # Storage degradation models
├── misc_errors/           # Additional error sources
├── sequencing/            # Sequencing platform simulators
└── read_processing/       # Read processing algorithms
```

### Adding New Methods

Create Python file in appropriate folder implementing the required interface. Pipeline discovers it automatically—specify method name in parameters to use immediately. Develop and test locally without core code modifications.

**Example:**
```python
# In dnabyte/encoding/my_encoder.py
class MyEncoder:
    def encode(self, data): return encoded_data, info
    def decode(self, data): return decoded_data, valid, info

# Usage:
params.encoding_method = "my_encoder"
```

### Plugin Discovery and Loading

`load_plugins.py` scans folders, imports modules, registers methods by name, and validates interfaces. This needs no core code modification

### Folder-Specific Guidelines

#### **`encoding/`** - Encoding Schemes
Must implement: `encode(data)`, `process(data)` and `decode(data)` methods

Examples: Reed-Solomon, Fountain codes, HEDGES, DNA Fountain

#### **`error_correction/`** - Error Correction Codes
Must implement: `encode(data)` and `decode(data)` with error correction

Examples: LDPC, Turbo codes, Polar codes

#### **`synthesis/`** - Synthesis Simulators
Must implement: `simulate(data)` returning data with synthesis errors

Examples: Empirical error models, platform-specific models (Twist, IDT, CustomArray)

#### **`storage/`** - Storage Degradation Models
Must implement: `simulate(data)` with storage condition parameters

Examples: Depurination models, strand break models, oxidative damage

#### **`sequencing/`** - Sequencing Platform Models
Must implement: `simulate(data)` returning sequenced reads

Examples: Illumina, Nanopore, PacBio HiFi models

#### **`misc_errors/`** - Additional Error Sources
Must implement: `simulate(data)` for other error types

Examples: PCR errors, contamination, sample loss

#### **`read_processing/`** - Read Processing Algorithms
Must implement: `process(data)` for quality filtering, clustering, consensus

Examples: Clustering algorithms, alignment methods, error correction

#### **`binarization/`** - Binarization Methods
Must implement: `binarize(data)` and `debinarize(data)`

Examples: Compression, encryption, format-specific handling

### Configuration Through Parameters

The `Params` object allows you to specify which implementations to use:

```python
params = Params()
params.encoding_method = "reed_solomon"
params.synthesis_method = "twist_bioscience"
params.storage_conditions = "room_temperature_1year"
params.sequencing_method = "illumina_novaseq"
params.error_methods = ["pcr_errors", "contamination"]
```


