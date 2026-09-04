#!/usr/bin/env python3
import os
from dnabyte.data_classes import Data
from dnabyte.binarize import Binarize
from dnabyte.params import Params

# Check file exists
file_path = './tests/testfiles/testfilesimsall.txt'
print(f"File exists: {os.path.exists(file_path)}")
print(f"File size: {os.path.getsize(file_path)} bytes")

# Load file
d = Data(file_paths=[file_path])
print(f"Data object attributes: {[a for a in dir(d) if not a.startswith('_')]}")

# Check binarize
params = Params(
    filename='testfilesimsall.txt',
    binarization_method='default',
    encoding_method='wukong',
    sequence_length=200,
)

bin = Binarize(params)
binary_code = bin.binarize(d)
print(f"After binarize - Data attributes: {[a for a in dir(binary_code) if not a.startswith('_')]}")
print(f"Binary code data attribute: {hasattr(binary_code, 'data')}")
if hasattr(binary_code, 'data'):
    print(f"Binary code data type: {type(binary_code.data)}")
    print(f"Binary code data length: {len(binary_code.data)}")
    if binary_code.data:
        print(f"First element type: {type(binary_code.data[0])}")
        print(f"First element length: {len(binary_code.data[0])}")
        print(f"First 100 chars: {binary_code.data[0][:100]}")
