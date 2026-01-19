#!/usr/bin/env python3
"""
Test suite for new InSilicoDNA class with bitarray packing.
"""

import sys
from bitarray import bitarray
from dnabyte.data_classes.new_insilicodna import InSilicoDNA

def test_basic_creation():
    """Test creating InSilicoDNA from DNA strings."""
    sequences = ['ATCG', 'GGCC', 'TATA']
    dna = InSilicoDNA(sequences)
    
    assert dna.num_sequences == 3
    assert dna.total_length == 12
    assert dna.average_length == 4.0
    print("✓ Basic creation works")

def test_packing_unpacking():
    """Test that packing and unpacking preserves data."""
    sequences = ['ATCGATCG', 'CCGGTTAA', 'AAAACCCC']
    dna = InSilicoDNA(sequences)
    
    # Verify unpacking returns original data
    assert dna.get_sequence(0) == 'ATCGATCG'
    assert dna.get_sequence(1) == 'CCGGTTAA'
    assert dna.get_sequence(2) == 'AAAACCCC'
    print("✓ Packing/unpacking preserves data")

def test_already_packed():
    """Test that already-packed data is handled correctly."""
    sequences = ['ATCG', 'GGCC']
    dna1 = InSilicoDNA(sequences)
    
    # Create from already-packed data
    dna2 = InSilicoDNA(dna1.data)
    
    assert dna2.get_sequence(0) == 'ATCG'
    assert dna2.get_sequence(1) == 'GGCC'
    print("✓ Already-packed data works (auto-detected)")

def test_indexing():
    """Test indexing and iteration."""
    sequences = ['ATCG', 'GGCC', 'TATA']
    dna = InSilicoDNA(sequences)
    
    # Test __getitem__
    assert dna[0] == 'ATCG'
    assert dna[1] == 'GGCC'
    assert dna[2] == 'TATA'
    
    # Test __len__
    assert len(dna) == 3
    
    # Test __iter__
    result = list(dna)
    assert result == sequences
    print("✓ Indexing and iteration work")

def test_add_remove():
    """Test adding and removing sequences."""
    sequences = ['ATCG', 'GGCC']
    dna = InSilicoDNA(sequences)
    
    # Add sequence
    dna.add_sequence('TATA')
    assert len(dna) == 3
    assert dna[2] == 'TATA'
    
    # Remove sequence
    dna.remove_sequence(1)
    assert len(dna) == 2
    assert dna[0] == 'ATCG'
    assert dna[1] == 'TATA'
    print("✓ Add/remove sequences work")

def test_validation():
    """Test validation of DNA sequences."""
    # Valid sequences
    valid = ['ATCG', 'GGCC']
    dna = InSilicoDNA(valid)
    assert dna.validate() == True
    
    # Invalid sequences
    try:
        invalid = ['ATCG', 'GGCX']  # X is invalid
        dna = InSilicoDNA(invalid)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "invalid nucleotides" in str(e)
    
    # Empty sequence
    try:
        empty = ['ATCG', '']
        dna = InSilicoDNA(empty)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "cannot be empty" in str(e)
    
    print("✓ Validation works")

def test_nucleotide_counts():
    """Test nucleotide counting."""
    sequences = ['ATCG', 'AAAA', 'TTTT']
    dna = InSilicoDNA(sequences)
    
    counts = dna.get_nucleotide_counts()
    assert counts['A'] == 5
    assert counts['T'] == 5
    assert counts['C'] == 1
    assert counts['G'] == 1
    print("✓ Nucleotide counts work")

def test_random():
    """Test random sequence generation."""
    dna = InSilicoDNA.random(10, 50)
    
    assert len(dna) == 10
    assert dna.num_sequences == 10
    assert dna.total_length == 500
    assert dna.average_length == 50.0
    
    # Verify all sequences are valid
    for seq in dna:
        assert len(seq) == 50
        assert all(c in 'ACGT' for c in seq)
    
    print("✓ Random generation works")

def test_complement():
    """Test complement and reverse complement."""
    # Test complement
    assert InSilicoDNA.complement('ATCG') == 'TAGC'
    assert InSilicoDNA.complement('AAAA') == 'TTTT'
    
    # Test reverse complement
    assert InSilicoDNA.reverse_complement('ATCG') == 'CGAT'
    assert InSilicoDNA.reverse_complement('AAAA') == 'TTTT'
    print("✓ Complement methods work")

def test_memory_efficiency():
    """Test memory efficiency of bitarray packing."""
    import sys
    
    # Create a large dataset
    n_sequences = 1000
    seq_length = 100
    
    # String-based storage
    string_sequences = [''.join(['ACGT'[i % 4] for i in range(seq_length)]) for _ in range(n_sequences)]
    string_size = sum(sys.getsizeof(seq) for seq in string_sequences)
    
    # Bitarray-based storage
    dna = InSilicoDNA(string_sequences)
    bitarray_size = sum(sys.getsizeof(seq) for seq in dna.data)
    
    savings = 100 * (1 - bitarray_size / string_size)
    print(f"✓ Memory efficiency: {savings:.1f}% savings (string: {string_size} bytes, bitarray: {bitarray_size} bytes)")
    
    # Should have significant savings (at least 20% due to Python overhead)
    assert savings > 20, f"Expected at least 20% savings, got {savings:.1f}%"

def run_all_tests():
    """Run all test functions."""
    test_functions = [
        test_basic_creation,
        test_packing_unpacking,
        test_already_packed,
        test_indexing,
        test_add_remove,
        test_validation,
        test_nucleotide_counts,
        test_random,
        test_complement,
        test_memory_efficiency,
    ]
    
    print("Running InSilicoDNA tests...")
    print("=" * 60)
    
    passed = 0
    failed = 0
    
    for test_func in test_functions:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"✗ {test_func.__name__} failed: {e}")
            failed += 1
    
    print("=" * 60)
    print(f"Tests passed: {passed}/{len(test_functions)}")
    
    if failed > 0:
        print(f"Tests failed: {failed}")
        sys.exit(1)
    else:
        print("All tests passed!")
        sys.exit(0)

if __name__ == '__main__':
    run_all_tests()
