#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test FASTA import functionality for NucleobaseCode
"""

import os
import tempfile
from dnabyte.data_classes.new_nucleobasecode import NucleobaseCode


def test_fasta_import_simple():
    """Test importing a simple FASTA file."""
    print("Test 1: Simple FASTA import")
    
    # Create temporary FASTA file
    fasta_content = """>sequence_1
ATCGATCG
>sequence_2
GCTAGCTA
>sequence_3
AAAATTTTCCCCGGGG
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_file = f.name
    
    try:
        # Import FASTA
        nbc = NucleobaseCode.from_fasta(fasta_file)
        
        # Verify
        assert nbc.num_codewords == 3, f"Expected 3 sequences, got {nbc.num_codewords}"
        
        decoded = nbc.to_string_format()
        assert decoded[0] == "ATCGATCG", f"Sequence 1 mismatch: {decoded[0]}"
        assert decoded[1] == "GCTAGCTA", f"Sequence 2 mismatch: {decoded[1]}"
        assert decoded[2] == "AAAATTTTCCCCGGGG", f"Sequence 3 mismatch: {decoded[2]}"
        
        print("✓ Simple FASTA import works")
        print(f"  Imported {nbc.num_codewords} sequences")
        print(f"  First sequence: {decoded[0]}")
        
    finally:
        os.unlink(fasta_file)


def test_fasta_import_multiline():
    """Test importing FASTA with multi-line sequences."""
    print("\nTest 2: Multi-line sequences")
    
    fasta_content = """>long_sequence
ATCGATCGATCG
GCTAGCTAGCTA
AAATTTCCCGGG
>another_sequence
TTTT
CCCC
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_file = f.name
    
    try:
        nbc = NucleobaseCode.from_fasta(fasta_file)
        decoded = nbc.to_string_format()
        
        expected_1 = "ATCGATCGATCGGCTAGCTAGCTAAAATTTCCCGGG"
        expected_2 = "TTTTCCCC"
        
        assert decoded[0] == expected_1, f"Multi-line sequence 1 mismatch"
        assert decoded[1] == expected_2, f"Multi-line sequence 2 mismatch"
        
        print("✓ Multi-line sequences handled correctly")
        print(f"  Sequence 1 length: {len(decoded[0])} bp")
        print(f"  Sequence 2 length: {len(decoded[1])} bp")
        
    finally:
        os.unlink(fasta_file)


def test_fasta_import_lowercase():
    """Test that lowercase sequences are converted to uppercase."""
    print("\nTest 3: Lowercase to uppercase conversion")
    
    fasta_content = """>lowercase_seq
atcgatcg
>mixed_case
AtCgAtCg
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_file = f.name
    
    try:
        nbc = NucleobaseCode.from_fasta(fasta_file)
        decoded = nbc.to_string_format()
        
        assert decoded[0] == "ATCGATCG", "Lowercase should be converted"
        assert decoded[1] == "ATCGATCG", "Mixed case should be converted"
        
        print("✓ Case conversion works correctly")
        
    finally:
        os.unlink(fasta_file)


def test_fasta_import_as_list():
    """Test importing with as_list=True option."""
    print("\nTest 4: Import as list structure")
    
    fasta_content = """>seq1
ATCG
>seq2
GCTA
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_file = f.name
    
    try:
        nbc = NucleobaseCode.from_fasta(fasta_file, as_list=True)
        decoded = nbc.to_string_format()
        
        # Should be list of lists
        assert isinstance(decoded[0], list), "Should be list structure"
        assert decoded[0][0] == "ATCG", "Nested sequence should match"
        assert decoded[1][0] == "GCTA", "Nested sequence should match"
        
        print("✓ List structure import works")
        print(f"  Structure: {type(decoded[0])}")
        
    finally:
        os.unlink(fasta_file)


def test_fasta_import_errors():
    """Test error handling."""
    print("\nTest 5: Error handling")
    
    # Test file not found
    try:
        NucleobaseCode.from_fasta("/nonexistent/file.fasta")
        assert False, "Should raise FileNotFoundError"
    except FileNotFoundError as e:
        print(f"✓ File not found error: {type(e).__name__}")
    
    # Test empty file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write("")
        empty_file = f.name
    
    try:
        try:
            NucleobaseCode.from_fasta(empty_file)
            assert False, "Should raise ValueError for empty file"
        except ValueError as e:
            print(f"✓ Empty file error: {type(e).__name__}")
    finally:
        os.unlink(empty_file)
    
    # Test invalid DNA characters
    fasta_content = """>invalid_seq
ATCGXYZ
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        invalid_file = f.name
    
    try:
        try:
            NucleobaseCode.from_fasta(invalid_file)
            assert False, "Should raise ValueError for invalid DNA"
        except ValueError as e:
            print(f"✓ Invalid DNA error: {type(e).__name__}")
    finally:
        os.unlink(invalid_file)


def test_fasta_encoding():
    """Verify that imported sequences are properly encoded."""
    print("\nTest 6: Bitarray encoding verification")
    
    fasta_content = """>test_seq
ATCG
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        fasta_file = f.name
    
    try:
        nbc = NucleobaseCode.from_fasta(fasta_file)
        
        # Check that internal data is bitarray
        from bitarray import bitarray
        assert isinstance(nbc.data[0], bitarray), "Should be encoded as bitarray"
        
        # Verify 2-bit encoding (4 nucleotides = 8 bits)
        assert len(nbc.data[0]) == 8, "ATCG should be 8 bits (2 bits per nucleotide)"
        
        print("✓ Sequences are properly encoded as bitarray")
        print(f"  ATCG encoded as: {nbc.data[0].to01()}")
        
    finally:
        os.unlink(fasta_file)


def run_all_tests():
    """Run all FASTA import tests."""
    print("=" * 60)
    print("Testing FASTA import for NucleobaseCode")
    print("=" * 60)
    
    try:
        test_fasta_import_simple()
        test_fasta_import_multiline()
        test_fasta_import_lowercase()
        test_fasta_import_as_list()
        test_fasta_import_errors()
        test_fasta_encoding()
        
        print("\n" + "=" * 60)
        print("✓ All FASTA import tests passed!")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


if __name__ == "__main__":
    import sys
    success = run_all_tests()
    sys.exit(0 if success else 1)
