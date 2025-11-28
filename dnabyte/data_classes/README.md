# DNAbyte Data Classes

This folder contains all data classes used in the DNAbyte pipeline for DNA-based data storage.

## Overview

The DNAbyte pipeline transforms data through three main stages, each represented by a specific data class. 

```
Data → BinaryCode → NucleobaseCode → InSilicoDNA → NucleobaseCode → BinaryCode → Data 
```

## Data Classes

### **Data** (`base.py`)
The foundation class for all data types in the DNAbyte pipeline.
- **Purpose**: Base class representing original file data
- **Usage**: `Data(file_paths)` or `Data.from_folder(folder_path)`
- **Features**: File validation, size calculation, path management

### **BinaryCode** (`binarycode.py`) 
Binary data representation consisting of a bitstream of '0' and '1' characters.
- **Purpose**: Digital representation of file data as binary sequences
- **Validation**: Ensures data contains only '0' and '1' characters
- **Usage**: `BinaryCode("101010110011")`
- **Features**: Bit manipulation, indexing, iteration capabilities

### **NucleobaseCode** (`nucleobasecode.py`)
Hierarchical codewords structure representing encoded nucleotide data for DNA synthesis.
- **Purpose**: Multi-level data structure where binary data has been encoded into DNA nucleotide sequences
- **Structure**: Nested lists with codewords containing nucleotide motifs
- **Usage**: `NucleobaseCode([[['A', 'C'], ['G', 'T']], [['A'], ['C', 'G']]])`
- **Features**: Hierarchical validation, codeword management, pooling structure representation

### **InSilicoDNA** (`insilicodna.py`)
Collection of synthesized DNA oligonucleotide sequences with realistic imperfections.
- **Purpose**: Represents physically synthesized DNA with synthesis errors and variations
- **Usage**: `InSilicoDNA(["ATCGATCG", "GCTAGCTA", "TTAATTAA"])`
- **Features**: DNA sequence validation, nucleotide counting, synthesis error modeling

## Pipeline Flow

1. **Binarization**: `Data` → `BinaryCode`
   - File content converted to binary representation
   - Handles compression and file format conversion

2. **Encoding**: `BinaryCode` → `NucleobaseCode` 
   - Binary data encoded into DNA nucleotide codewords
   - Creates hierarchical structure for synthesis pooling

3. **Synthesis Simulation**: `NucleobaseCode` → `InSilicoDNA`
   - Simulates physical DNA synthesis process
   - Introduces realistic synthesis errors and variations

4. **Storage Simulation**: `InSilicoDNA`→ `InSilicoDNA`
   - Simulates physical DNA storage process
   - Introduces strand breaks based on probabilistic model

5. **Sequencing Simulation**: `InSilicoDNA`→ `InSilicoDNA`
   - Simulates physical sequencing process
   - Introduces different errors based on probabilistic model

6. **Processing**: `InSilicoDNA` → `NucleobaseCode`
   - Error correction and sequence reconstruction
   - Recovers original codeword structure

7. **Decoding**: `NucleobaseCode` → `BinaryCode`
   - Converts nucleotide codewords back to binary data
   - Reverses the encoding process

8. **Restoration**: `BinaryCode` → `Data`
   - Reconstructs original files from binary data
   - Handles decompression and file format restoration


## Usage Examples

### Direct Instantiation
```python
from dnabyte.dataclasses import Data, BinaryCode, NucleobaseCode, InSilicoDNA

# Create data from files
data = Data(['/path/to/file1.txt', '/path/to/file2.txt'])

# Create binary data
binary_code = BinaryCode("101010110011")

# Create hierarchical encoded data
nucleobase_code = NucleobaseCode([
    [['A', 'T'], ['G', 'C']],  # First codeword
    [['C', 'G'], ['A', 'T']]   # Second codeword
])

# Create synthesized DNA sequences
dna = InSilicoDNA(["ATCGATCG", "GCTAGCTA", "TTAATTAA"])
```
