# DNAbyte
**A comprehensive Python toolkit for DNA data storage simulation and analysis**

DNAbyte provides a complete framework for simulating the entire end-to-end pipeline of DNA-based data storage, from raw files to DNA sequences and back. Whether you're researching novel encoding schemes, testing error correction methods, or analyzing storage durability.

---

## Pipeline Overview

DNAbyte covers the complete DNA storage workflow:

1. **Binarization**: Compress and convert input files into binary sequences
2. **Encoding**: Transform binary data into DNA sequences (including error correction)
3. **Synthesis**: Simulate DNA synthesis with configurable error models (MESA, assembly-based)
4. **Storage**: Model degradation over time under various environmental conditions
5. **Sequencing**: Simulate sequencing technologies (Illumina, Nanopore, IID)
6. **Processing**: Align and correct sequenced reads
7. **Decoding**: Recover original binary data with error correction
8. **Restoration**: Decompress and restore original files

![Overview](figures/overview.png)

---

## Documentation Structure

* **[Usage Guide](usage.md)** - Complete examples and parameter documentation
* **[Simulation Guide](simulation.md)** - Step-by-step simulation tutorials  
* **[Core Concepts](concepts.md)** - Understanding encoding techniques and algorithms
* **[API Reference](reference.md)** - Comprehensive class and function documentation


---
## Installation

```bash
# Clone the repository
git clone https://github.com/Eko-Refugium/DNAbyte.git
cd DNAbyte

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

---

## Quick Start

```python
from dnabyte.params import Params
from dnabyte.data_classes.base import Data
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode

# Set up parameters
params = Params(
    file_paths=['data.txt'],
    encoding_method='max_density',
    binarization_method='compressed',
    outer_error_correction='reedsolomon',
    reed_solo_percentage=0.8
)

# Run the pipeline
data = Data(file_paths=params.file_paths)
binarizer = Binarize(params)
encoder = Encode(params)

binary_code = binarizer.binarize(data)
nucleobase_code, info = encoder.encode(binary_code)
```

---



## Citation

If you use DNAbyte in your research, please cite:

```bibtex
@software{dnabyte2025,
  title = {DNAbyte: A Framework for End-to-end DNA Data Storage Simulation},
  author = {Kaya Wernhart and Fabian Schroeder},
  year = {2025},
  url = {https://github.com/Eko-Refugium/DNAbyte}
}
```

---

## Acknowledgments

Eko Refugium is part of the **MI-DNA DISC consortium** ([www.midnadisc.eu](https://www.midnadisc.eu)).

MI-DNA DISC was funded by the EU Commission in the framework of the **Horizon Europe – EIC Pathfinder Challenges** programme.

**Grant Agreement**: 101115215

![EU Logo](assets/EU_logo.png)

---

## License

See [LICENSE](../LICENSE) for details.