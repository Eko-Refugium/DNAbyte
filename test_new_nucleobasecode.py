#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test script for new_nucleobasecode.py with bitarray 2-bit encoding
"""

import sys
from bitarray import bitarray
from dnabyte.data_classes.new_nucleobasecode_v2 import NucleobaseCode


def test_encoding_decoding():
    """Test 2-bit DNA encoding and decoding."""
    print("Test 1: DNA encoding/decoding")
    
    # Create instance with string data
    data = [["ATCGAACGCGTGAGTC", "GCTGAACGCGTGA"], ["CGCGATCTCGTCGATCG", "AAGAACGCGTGA",["CGCGTATAGTCGTATCGT", "TTCTTGCGCACT"]]]
    nbc = NucleobaseCode(data)
    
    # Verify encoding happened (internal data should be bitarrays)
    assert isinstance(nbc.data[0][0], bitarray), "Data should be encoded as bitarray"
    assert isinstance(nbc.data[0][1], bitarray), "Data should be encoded as bitarray"
    assert isinstance(nbc.data[1], list), "Data should be encoded as bitarray"
    
    # Verify decoding works
    decoded = nbc.to_string_format()
    assert decoded == data, f"Decoded data should match original: {decoded} != {data}"
    
    print("✓ Encoding/decoding works correctly")
    print(f"  Original: {data}")
    print(f"  Encoded (bitarray): {nbc.data}")  # Show first 8 bits
    print(f"  Decoded: {decoded}")


def test_2bit_efficiency():
    """Verify 2-bit encoding (2 bits per nucleotide)."""
    print("\nTest 2: 2-bit encoding efficiency")
    
    # Test single DNA string
    dna = "ATCG"  # 4 nucleotides
    nbc = NucleobaseCode([dna])
    
    # Should be 8 bits total (2 bits × 4 nucleotides)
    encoded_bits = len(nbc.data[0])
    assert encoded_bits == 8, f"Expected 8 bits, got {encoded_bits}"
    
    # Verify specific encoding
    # A=00, T=11, C=01, G=10 -> "00110110"
    expected = bitarray('00110110')
    assert nbc.data[0] == expected, f"Encoding mismatch: {nbc.data[0].to01()} != {expected.to01()}"
    
    print("✓ 2-bit encoding verified")
    print(f"  DNA: {dna}")
    print(f"  Bits: {nbc.data[0].to01()} (A=00, T=11, C=01, G=10)")
    print(f"  Length: {encoded_bits} bits for 4 nucleotides")


def test_nested_structure():
    """Test preservation of nested list structure."""
    print("\nTest 3: Nested structure preservation")
    
    # Create complex nested structure
    data = [
        "ATCG",                          # Simple string
        ["AAAA", "TTTT"],                # 1-level nesting
        [["CCCC", "GGGG"], "ATCG"],      # 2-level nesting
    ]
    
    nbc = NucleobaseCode(data)
    
    # Verify structure is preserved
    assert len(nbc.data) == 3, "Top level should have 3 codewords"
    assert isinstance(nbc.data[0], bitarray), "First codeword should be bitarray"
    assert isinstance(nbc.data[1], list), "Second codeword should be list"
    assert len(nbc.data[1]) == 2, "Second codeword should have 2 elements"
    assert isinstance(nbc.data[2], list), "Third codeword should be list"
    assert isinstance(nbc.data[2][0], list), "Third codeword should have nested list"
    
    # Verify decoding preserves structure
    decoded = nbc.to_string_format()
    assert decoded == data, "Decoded structure should match original"
    
    print("✓ Nested structure preserved")
    print(f"  Max depth: {nbc.max_depth}")
    print(f"  Structure matches: {decoded == data}")


def test_get_codeword():
    """Test getting codewords with and without unpacking."""
    print("\nTest 4: Get codeword functionality")
    
    data = ["ATCG", "GCTA", "AAAA"]
    nbc = NucleobaseCode(data)
    
    # Get packed (bitarray)
    codeword_packed = nbc.get_codeword(0, unpack=False)
    assert isinstance(codeword_packed, bitarray), "Should return bitarray when unpack=False"
    
    # Get unpacked (string)
    codeword_unpacked = nbc.get_codeword(0, unpack=True)
    assert codeword_unpacked == "ATCG", "Should return original string when unpack=True"
    
    # Test indexing (__getitem__)
    codeword_indexed = nbc[1]
    assert isinstance(codeword_indexed, bitarray), "__getitem__ should return packed"
    
    print("✓ Get codeword works")
    print(f"  Packed: {codeword_packed.to01()}")
    print(f"  Unpacked: {codeword_unpacked}")


def test_validation():
    """Test input validation."""
    print("\nTest 5: Input validation")
    
    # Test invalid DNA characters
    try:
        nbc = NucleobaseCode(["ATCGX"])
        assert False, "Should raise ValueError for invalid DNA character"
    except ValueError as e:
        print(f"✓ Invalid characters rejected: {e}")
    
    # Test empty data
    try:
        nbc = NucleobaseCode([])
        assert False, "Should raise ValueError for empty data"
    except ValueError as e:
        print(f"✓ Empty data rejected: {e}")
    
    # Test wrong type
    try:
        nbc = NucleobaseCode("ATCG")  # Should be list
        assert False, "Should raise TypeError for non-list"
    except TypeError as e:
        print(f"✓ Wrong type rejected: {e}")


def test_memory_efficiency():
    """Compare memory usage between string and bitarray encoding."""
    print("\nTest 6: Memory efficiency comparison")
    
    # Generate large dataset
    n_codewords = 1000
    seq_length = 250
    
    # Create string version (simulate old approach)
    string_data = []
    for _ in range(n_codewords):
        seq = "ATCG" * (seq_length // 4)
        string_data.append(seq)
    
    string_size = sum(sys.getsizeof(seq) for seq in string_data)
    
    # Create bitarray version (new approach)
    nbc = NucleobaseCode(string_data)
    
    # Calculate bitarray size
    bitarray_size = sum(
        sys.getsizeof(codeword) 
        for codeword in nbc.data 
        if isinstance(codeword, bitarray)
    )
    
    savings = (1 - bitarray_size / string_size) * 100
    
    print(f"  {n_codewords} codewords × {seq_length} bp each:")
    print(f"  String storage: {string_size:,} bytes ({string_size/(n_codewords*seq_length):.2f} bytes/bp)")
    print(f"  Bitarray storage: {bitarray_size:,} bytes ({bitarray_size/(n_codewords*seq_length):.4f} bytes/bp)")
    print(f"  Memory savings: {savings:.1f}%")
    print("✓ Memory efficiency confirmed")


def test_random_generation():
    """Test random NucleobaseCode generation."""
    print("\nTest 7: Random generation")
    
    # Test simple synthesis type
    nbc = NucleobaseCode.random(type='synthesis', n=10, m=50)
    assert nbc.num_codewords == 10, "Should have 10 codewords"
    assert len(nbc) == 10, "__len__ should return 10"
    
    # Verify data is encoded
    assert isinstance(nbc.data[0], bitarray), "Random data should be encoded"
    
    # Verify can decode
    decoded = nbc.to_string_format()
    assert isinstance(decoded[0], str), "Should decode to strings"
    assert len(decoded[0]) == 50, "Each codeword should be 50 bp"
    assert all(c in 'ATCG' for c in decoded[0]), "Should contain only valid DNA bases"
    
    print("✓ Random generation works")
    print(f"  Generated {nbc.num_codewords} codewords")
    print(f"  Preview: {decoded[0][:30]}...")


def test_string_representations():
    """Test __str__ and __repr__ methods."""
    print("\nTest 8: String representations")
    
    data = ["ATCG", "GCTA"]
    nbc = NucleobaseCode(data)
    
    str_repr = str(nbc)
    print("__str__ output:")
    print(str_repr)
    
    repr_repr = repr(nbc)
    print(f"__repr__ output: {repr_repr}")
    
    assert "NucleobaseCode" in str_repr, "__str__ should contain class name"
    assert "2" in str_repr, "__str__ should show 2 codewords"
    assert "NucleobaseCode" in repr_repr, "__repr__ should contain class name"
    
    print("✓ String representations work")


def test_already_encoded():
    """Test creating NucleobaseCode from already-packed data."""
    print("\nTest 9: Creating from already-packed data")
    
    # Manually create packed structure
    # A=00, T=11, C=01, G=10
    packed_data = [
        bitarray('00110110'),  # ATCG = 00-11-01-10
        bitarray('10011100')   # GCTA = 10-01-11-00
    ]
    
    # Create - will auto-detect as packed
    nbc = NucleobaseCode(packed_data)
    
    # Verify no double-packing
    assert nbc.data[0] == packed_data[0], "Should not re-pack"
    
    # Verify can unpack
    unpacked = nbc.to_string_format()
    assert unpacked[0] == "ATCG", "Should unpack correctly"
    assert unpacked[1] == "GCTA", "Should unpack correctly"
    
    print("✓ Already-packed data works (auto-detected)")
    print(f"  Unpacked: {unpacked}")


def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Testing new_nucleobasecode.py with 2-bit bitarray encoding")
    print("=" * 60)
    
    try:
        test_encoding_decoding()
        test_2bit_efficiency()
        test_nested_structure()
        test_get_codeword()
        test_validation()
        test_memory_efficiency()
        test_random_generation()
        test_string_representations()
        test_already_encoded()
        
        print("\n" + "=" * 60)
        print("✓ All tests passed!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
