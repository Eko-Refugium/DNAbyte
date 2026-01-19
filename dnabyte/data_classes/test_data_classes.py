#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test script for new_binarycode.py with bitarray implementation
"""

import sys
import os
from bitarray import bitarray
from dnabyte.library import Library
from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.new_binarycode import BinaryCode as NewBinaryCode
from dnabyte.data_classes.new_nucleobasecode_v2 import NucleobaseCode
from dnabyte.data_classes.new_insilicodna import InSilicoDNA as NewInSilicoDNA
from dnabyte.data_classes.insilicodna import InSilicoDNA


from dnabyte.params import Params
from dnabyte.binarize import Binarize  


# ## IDENTIFY PROBLEM WITH BINARYCODE SIZE ##

# # Print file size
# file_size_1 = os.path.getsize('./tests/testfiles/Bohemian_Rhapsody_Lyrics.txt')
# file_size_2 = os.path.getsize('./tests/testfiles/Feelin_Good_Lyrics.txt')
# print(f"Size: {file_size_1+file_size_2} bytes ({(file_size_1 + file_size_2) / 1024:.2f} KB)")

# file_paths = ['./tests/testfiles/Bohemian_Rhapsody_Lyrics.txt', './tests/testfiles/Feelin_Good_Lyrics.txt']
# params = Params(binarization_method='compressed', file_paths=file_paths)

# data = Data(file_paths=file_paths)
# bin = Binarize(params=params)
# binary_code = bin.binarize(data)

# # Print binary_code size
# binary_code_size = sys.getsizeof(binary_code.data)
# print(f"Binary code size: {binary_code_size} bytes ({binary_code_size / 1024:.2f} KB)")
# print(f"Compression ratio: {(file_size_1 + file_size_2) / binary_code_size:.2f}x")


# ## 

# random_binary_code = BinaryCode.random(n=1000)
# #print(random_binary_code.data)

# random_binary_code_size = sys.getsizeof(random_binary_code.data)
# print(f"Binary code size: {random_binary_code_size} bytes ({random_binary_code_size / 1024:.2f} KB)")

# random_binary_code_new = NewBinaryCode.random(n=1000)
# print(random_binary_code_new.data)
# random_binary_code_new_size = sys.getsizeof(random_binary_code_new.data)
# print(f"New Binary code size: {random_binary_code_new_size} bytes ({random_binary_code_new_size / 1024:.2f} KB)")


# ## TEST NUCLEOBASECODE V2 PACKING ##

# nbc = [['ATCG', 'GCTA', 'TTAA'], [[['CCGG', 'AATT'], ['GGCC', 'TTAA']], 'CGTA']]
# nbc = NucleobaseCode(nbc)

# # packed
# print(f"  _structure: {nbc._structure}")
# print(f"  _packed_data: {nbc._packed_data.to01()}")

# # unpacked
# print(nbc.unpack())

# data_lib_indices = [[34, 197, 65], [[[113, 254], [77, 171]], 42]]
# lib = Library(structure='linear_assembly', filename='./tests/testlibraries/20bp_Lib.csv')
# nbc = NucleobaseCode(data_lib_indices, library=lib)


# # packed
# print(f"  _structure: {nbc._structure}")
# print(f"  _packed_data: {nbc._packed_data.to01()}")

# # unpacked
# print(nbc.unpack())
# print(lib.library[34])


# ## TEST INSILICODNA PACKING ##

is_dna = ['ATCG', 'GCTA', 'TTAA', 'CCGG', 'AATT', 'GGCC', 'TTAA', 'CGTA']
isdna_new = NewInSilicoDNA(is_dna)
isdna_old = InSilicoDNA(is_dna)

print("Old InSilicoDNA data:")
print(isdna_old.data) 

print("\nNew InSilicoDNA data:")
print(isdna_new._packed_data)
print(isdna_new._sequence_lengths)



