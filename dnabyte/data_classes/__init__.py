"""
Data classes for DNAbyte - representations of data at different stages of encoding/decoding.
"""

from .binarycode import BinaryCode
from .nucleobasecode import NucleobaseCode
from .insilicodna import InSilicoDNA

__all__ = [
    "BinaryCode",
    "NucleobaseCode", 
    "InSilicoDNA",
]
