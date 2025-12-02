"""
Data classes for DNAbyte - representations of data at different stages of encoding/decoding.
"""

from .base import Data
from .binarycode import BinaryCode
from .nucleobasecode import NucleobaseCode
from .insilicodna import InSilicoDNA

__all__ = [
    "Data",
    "BinaryCode",
    "NucleobaseCode", 
    "InSilicoDNA",
]
