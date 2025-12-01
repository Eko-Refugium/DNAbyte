# DNAbyte
A comprehensive Python toolkit for DNA-based data storage.

DNAbyte provides an end-to-end simulation framework for DNA data storage systems. The toolkit comprises five core components: (i) binarization, (ii) encoding, (iii) synthesis/storage/sequencing simulation, (iv) decoding, and (v) file restoration. Each component can be used independently or combined into a seamless pipeline, enabling researchers to simulate and benchmark different DNA data storage methodologies.

DNAbyte was developed under the EIC Pathfinder Challenge MI-DNA DISC (www.midnadisc.eu). This project aims to bring a low-cost, energy-efficient, and fast data driven approach that can write, edit, store, and retrieve DNA-based data, more efficiently compared to current technologies. It was funded by the EU Commission in the framework of the Horizon Europe – EIC Pathfinder Challenges programme under the Grant Agreement 101115215.

## Installation

DNAbyte uses [uv](https://github.com/astral-sh/uv), a fast Python package and project manager.

### Install uv

First, install uv if you don't have it:

```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Or with Homebrew
brew install uv
```

### Create a new project and install DNAbyte

```bash
# Create a new project
uv init my-dna-project
cd my-dna-project

# Add DNAbyte from GitHub
uv add "dnabyte @ git+ssh://git@github.com/Eko-Refugium/DNAbyte.git"

# Sync dependencies
uv sync

# Run your Python scripts
uv run python your_script.py
```

### Add DNAbyte to an existing project

```bash
cd your-existing-project
uv add "dnabyte @ git+ssh://git@github.com/Eko-Refugium/DNAbyte.git"
uv sync
```

### Verification

Verify the installation by importing DNAbyte:

```python
import dnabyte
from dnabyte.data_classes.base import Data

print("DNAbyte installed successfully!")
```

### Important Notes

- **Private Submodules**: DNAbyte uses private Git submodules for some encoding algorithms. Access to the repository and its submodules is required.
- **Dependencies**: All required dependencies (scipy, numpy, biopython, etc.) will be installed automatically.`


