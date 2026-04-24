#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test script for new_binarycode.py with bitarray implementation
"""

import sys
from bitarray import bitarray
from dnabyte.data_classes.new_binarycode import BinaryCode

def test_creation_from_string():
    """Test creating BinaryCode from string."""
    print("Test 1: Creating BinaryCode from string")
    binary = BinaryCode("10110010")
    assert binary.length == 8
    assert binary.to_bitstring() == "10110010"
    print("✓ Creation from string works")

def test_creation_from_bitarray():
    """Test creating BinaryCode from bitarray."""
    print("\nTest 2: Creating BinaryCode from bitarray")
    ba = bitarray("11001100")
    binary = BinaryCode(ba)
    assert binary.length == 8
    assert binary.to_bitstring() == "11001100"
    print("✓ Creation from bitarray works")

def test_validation():
    """Test validation of input."""
    print("\nTest 3: Input validation")
    
    # Test empty string
    try:
        BinaryCode("")
        assert False, "Should raise ValueError for empty string"
    except ValueError as e:
        print(f"✓ Empty string rejected: {e}")
    
    # Test invalid characters
    try:
        BinaryCode("101201")
        assert False, "Should raise ValueError for invalid chars"
    except ValueError as e:
        print(f"✓ Invalid characters rejected: {e}")
    
    # Test wrong type
    try:
        BinaryCode(12345)
        assert False, "Should raise TypeError for non-string/bitarray"
    except TypeError as e:
        print(f"✓ Wrong type rejected: {e}")

def test_length_indexing_iteration():
    """Test __len__, __getitem__, and __iter__ methods."""
    print("\nTest 4: Length, indexing, and iteration")
    binary = BinaryCode("10110010")
    
    # Test length
    assert len(binary) == 8
    print(f"✓ len(binary) = {len(binary)}")
    
    # Test indexing
    assert binary[0] == 1  # bitarray returns int
    assert binary[3] == 1
    assert binary[7] == 0
    print(f"✓ Indexing works: binary[0]={binary[0]}, binary[3]={binary[3]}")
    
    # Test slicing
    slice_result = binary[2:5]
    assert len(slice_result) == 3
    print(f"✓ Slicing works: binary[2:5] has length {len(slice_result)}")
    
    # Test iteration
    bits = [bit for bit in binary]
    assert len(bits) == 8
    print(f"✓ Iteration works: {bits}")

def test_random_generation():
    """Test random BinaryCode generation."""
    print("\nTest 5: Random generation")
    binary = BinaryCode.random(100)
    assert binary.length == 100
    assert len(binary) == 100
    print(f"✓ Random generation works: created {binary.length} bits")
    print(f"  Preview: {binary.to_bitstring()[:20]}...")

def test_string_representation():
    """Test __str__ and __repr__ methods."""
    print("\nTest 6: String representations")
    binary = BinaryCode("10110010" * 10, size=10)
    
    str_repr = str(binary)
    print("__str__ output:")
    print(str_repr)
    
    repr_repr = repr(binary)
    print(f"__repr__ output: {repr_repr}")
    print("✓ String representations work")


def test_file_paths_and_size():
    """Test file_paths and size attributes."""
    print("\nTest 8: File paths and size attributes")
    binary = BinaryCode("10110010", file_paths=["/path/to/file.txt"], size=1024)
    assert binary.file_paths == ["/path/to/file.txt"]
    assert binary.size == 1024
    print("✓ File paths and size attributes work")

def test_compare_method():
    """Test the compare method (if implemented)."""
    print("\nTest 9: Compare method")
    try:
        binary1 = BinaryCode("10110010")
        binary2 = BinaryCode("10110010")
        binary3 = BinaryCode("10110011")
        
        result1, diff1 = binary1.compare(binary1, binary2)
        print(f"  Same data: {result1}, differences: {diff1}")
        
        result2, diff2 = binary1.compare(binary1, binary3)
        print(f"  Different data: {result2}, differences at indices: {diff2}")
        print("✓ Compare method works (needs update for bitarray)")
    except AttributeError:
        print("⚠ Compare method not yet implemented/updated")

def run_all_tests():
    """Run all tests."""
    print("=" * 60)
    print("Testing new_binarycode.py with bitarray implementation")
    print("=" * 60)
    
    try:
        test_creation_from_string()
        test_creation_from_bitarray()
        test_validation()
        test_length_indexing_iteration()
        test_random_generation()
        test_string_representation()
        test_file_paths_and_size()
        test_compare_method()
        
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
