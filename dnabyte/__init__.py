"""
DNAbyte - A toolkit for encoding and decoding data for DNA-based data storage
"""

__version__ = "0.9.0"

# Import main modules for easy access
from . import binarize
from . import encode
from . import store
from . import synthesize
from . import params
from . import oligo
from . import oligopool
from . import sequence
from . import library

# Import submodules
from . import binarization
from . import data_classes
from . import encoding
from . import error_correction
from . import sequencing
from . import storage
from . import synthesis

# Import commonly used data classes for convenience
from .data_classes import BinaryCode, NucleobaseCode, InSilicoDNA

__all__ = [
    "binarize",
    "encode", 
    "store",
    "synthesize",
    "params",
    "oligo",
    "oligopool",
    "sequence",
    "library",
    "binarization",
    "data_classes",
    "encoding",
    "error_correction",
    "sequencing",
    "storage",
    "synthesis",
    "BinaryCode",
    "NucleobaseCode",
    "InSilicoDNA",
]
