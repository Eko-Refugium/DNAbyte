#!/usr/bin/env python3
"""
Test suite for new_nucleobasecode_v2.py with single-bitarray storage.
"""

import sys
from bitarray import bitarray
from dnabyte.data_classes.new_nucleobasecode_v2 import NucleobaseCode


def test_basic_init():
    """Test creating NucleobaseCode from nested DNA strings."""
    print("Test 1: Basic creation and unpacking")
    
    data = [
        ['ATCGATCG', 'GGCCGGCC'],
        ['TTAATTAA', 'CCGGCCGG', 'AACGCGCAACCCC']
    ]
    
    nbc = NucleobaseCode(data)
    
    # Top level has 2 sublists
    assert nbc.num_codewords == 2
    assert nbc.max_depth > 0
    
    # Unpack and verify
    unpacked = nbc.unpack()
    assert unpacked == data, f"Unpacking failed: {unpacked} != {data}"
    
    print("✓ Basic creation and unpacking works")
    print(f"  Codewords: {nbc.num_codewords}")
    print(f"  Max depth: {nbc.max_depth}")
    print(f"  Packed size: {len(nbc._packed_data)} bits")
    print(f"  _structure: {nbc._structure}")
    print(f"  _packed_data: {nbc._packed_data.to01()}")
    print(f"  _structure: {nbc._structure}")
    print(f"  _packed_data: {nbc._packed_data.to01()}")



def test_library_indices_init():
    """Test library index packing mode."""
    print("\nTest 12: Library index packing")
    
    from dnabyte.library import Library
    import os
    
    # Skip if library file doesn't exist
    lib_path = './tests/testlibraries/20bp_Lib.csv'
    if not os.path.exists(lib_path):
        print("⊘ Library index test skipped (20bp_Lib.csv not found)")
        return
    
    # Load library
    lib = Library(structure='linear_assembly', filename=lib_path)
    
    # Create data with indices
    data = [
        [[0, 1, 2], [3, 4, 5]],
        [[6, 7, 8]]
    ]
    
    nbc = NucleobaseCode.from_library_indices(data, lib)
    
    assert nbc.pack_mode == 'library_indices'
    assert nbc.num_codewords == 2
    
    # Unpack should give library sequences
    unpacked = nbc.unpack()
    assert isinstance(unpacked[0][0][0], str)
    assert len(unpacked[0][0]) == 3
    
    print("✓ Library index packing works")
    print(f"  Index bits: {nbc.index_bits}")
    print(f"  Packed size: {len(nbc._packed_data)} bits")
    print(f"  _structure: {nbc._structure}")
    print(f"  _packed_data: {nbc._packed_data.to01()}")







def test_single_bitarray_storage():
    """Test that all data is stored in a single bitarray."""
    print("\nTest 2: Single bitarray storage")
    
    data = [
        ['ATCG', 'GGCC'],
        ['TTAA', 'CCGG']
    ]
    
    nbc = NucleobaseCode(data)
    
    # Verify single bitarray
    assert isinstance(nbc._packed_data, bitarray), "Should use single bitarray"
    assert isinstance(nbc._structure, list), "Should have structure metadata"
    
    # Total should be 16 nucleotides * 2 bits = 32 bits
    expected_bits = 16 * 2
    assert len(nbc._packed_data) == expected_bits, f"Expected {expected_bits} bits, got {len(nbc._packed_data)}"
    
    print("✓ Single bitarray storage confirmed")
    print(f"  Single bitarray: {len(nbc._packed_data)} bits")
    print(f"  Structure tokens: {len(nbc._structure)}")


def test_nested_structure_preservation():
    """Test that complex nested structures are preserved."""
    print("\nTest 3: Nested structure preservation")
    
    data = [
        'ATCG',
        ['AAAA', 'TTTT'],
        [['CCCC', 'GGGG'], 'ATCG'],
    ]
    
    nbc = NucleobaseCode(data)
    unpacked = nbc.unpack()
    
    assert unpacked == data, "Structure not preserved"
    assert len(unpacked) == 3
    assert isinstance(unpacked[0], str)
    assert isinstance(unpacked[1], list)
    assert isinstance(unpacked[2], list)
    assert isinstance(unpacked[2][0], list)
    
    print("✓ Nested structure preserved correctly")
    print(f"  Max depth: {nbc.max_depth}")


def test_2bit_encoding():
    """Verify 2-bit DNA encoding."""
    print("\nTest 4: 2-bit DNA encoding")
    
    # A=00, T=11, C=01, G=10
    dna = "ATCG"
    nbc = NucleobaseCode([dna])
    
    # Should be 8 bits (4 nucleotides * 2 bits)
    assert len(nbc._packed_data) == 8
    
    # Verify encoding: ATCG = 00-11-01-10
    expected_bits = '00110110'
    assert nbc._packed_data.to01() == expected_bits
    
    # Verify unpacking
    unpacked = nbc.unpack()
    assert unpacked[0] == dna
    
    print("✓ 2-bit encoding verified")
    print(f"  DNA: {dna}")
    print(f"  Bits: {nbc._packed_data.to01()}")


def test_indexing():
    """Test __getitem__ and __len__."""
    print("\nTest 5: Indexing operations")
    
    data = [['ATCG', 'GGCC'], ['TTAA', 'CCGG'], ['AAAA']]
    nbc = NucleobaseCode(data)
    
    # Test __len__ - should return 3 top-level codewords
    assert len(nbc) == 3
    
    # Test __getitem__
    assert nbc[0] == ['ATCG', 'GGCC']
    assert nbc[1] == ['TTAA', 'CCGG']
    assert nbc[2] == ['AAAA']
    
    print("✓ Indexing works")


def test_iteration():
    """Test __iter__."""
    print("\nTest 6: Iteration")
    
    data = [['ATCG'], ['GGCC'], ['TTAA']]
    nbc = NucleobaseCode(data)
    
    result = list(nbc)
    assert result == data
    
    print("✓ Iteration works")


def test_validation():
    """Test data validation."""
    print("\nTest 7: Data validation")
    
    # Valid data
    valid = [['ATCG', 'GGCC']]
    nbc = NucleobaseCode(valid)
    assert nbc.validate() == True
    
    # Invalid nucleotide
    try:
        invalid = [['ATCGX']]
        nbc = NucleobaseCode(invalid)
        assert False, "Should raise ValueError"
    except ValueError as e:
        assert "invalid nucleotides" in str(e).lower()
    
    # Empty list
    try:
        empty = [[]]
        nbc = NucleobaseCode(empty)
        assert False, "Should raise ValueError"
    except ValueError:
        pass
    
    # Wrong type
    try:
        wrong = "ATCG"
        nbc = NucleobaseCode(wrong)
        assert False, "Should raise TypeError"
    except TypeError:
        pass
    
    print("✓ Validation works")


def test_already_packed():
    """Test creating from already-packed format."""
    print("\nTest 8: Creating from already-packed data")
    
    data = [['ATCG', 'GGCC']]
    nbc1 = NucleobaseCode(data)
    
    # Create from packed format
    nbc2 = NucleobaseCode(nbc1.data)
    
    unpacked = nbc2.unpack()
    assert unpacked == data
    
    print("✓ Already-packed data works")


def test_memory_efficiency():
    """Compare memory with multiple-bitarray approach."""
    print("\nTest 9: Memory efficiency")
    
    # Create nested structure with many small sequences
    n_codewords = 100
    pools_per_cw = 3
    seqs_per_pool = 2
    seq_length = 30  # Short sequences where overhead matters
    
    data = []
    for _ in range(n_codewords):
        codeword = []
        for _ in range(pools_per_cw):
            pool = []
            for _ in range(seqs_per_pool):
                seq = 'ACGT' * (seq_length // 4)
                pool.append(seq)
            codeword.append(pool)
        data.append(codeword)
    
    nbc = NucleobaseCode(data)
    
    # Calculate actual storage
    packed_size = sys.getsizeof(nbc._packed_data)
    structure_size = sys.getsizeof(nbc._structure) + sum(
        sys.getsizeof(item) for item in nbc._structure if isinstance(item, int)
    )
    total_size = packed_size + structure_size
    
    # Calculate string storage
    def calculate_string_size(d):
        if isinstance(d, str):
            return sys.getsizeof(d)
        elif isinstance(d, list):
            return sys.getsizeof(d) + sum(calculate_string_size(item) for item in d)
        return 0
    
    string_size = calculate_string_size(data)
    
    savings = 100 * (1 - total_size / string_size)
    
    print(f"  Data: {n_codewords} codewords, {seq_length}bp sequences")
    print(f"  String storage: {string_size:,} bytes")
    print(f"  New storage:    {total_size:,} bytes")
    print(f"    - Packed:     {packed_size:,} bytes (single bitarray)")
    print(f"    - Structure:  {structure_size:,} bytes")
    print(f"  Memory savings: {savings:.1f}%")
    
    assert savings > 50, f"Expected >50% savings, got {savings:.1f}%"
    print("✓ Memory efficiency confirmed")


def test_complement():
    """Test complement and reverse_complement."""
    print("\nTest 10: Complement methods")
    
    assert NucleobaseCode.complement('ATCG') == 'TAGC'
    assert NucleobaseCode.complement('AAAA') == 'TTTT'
    assert NucleobaseCode.reverse_complement('ATCG') == 'CGAT'
    
    print("✓ Complement methods work")


def test_random_generation():
    """Test random data generation."""
    print("\nTest 11: Random generation")
    
    # Simple structure: list with 3 sublists, each having 2 sequences
    nbc = NucleobaseCode.random([2, 2, 2], min_length=20, max_length=30)
    
    # Should have 3 top-level codewords
    assert nbc.num_codewords == 3
    unpacked = nbc.unpack()
    assert len(unpacked) == 3
    assert all(len(cw) == 2 for cw in unpacked)
    
    # Verify sequences are valid DNA
    for cw in unpacked:
        for seq in cw:
            assert isinstance(seq, str)
            assert 20 <= len(seq) <= 30
            assert all(c in 'ACGT' for c in seq)
    
    print("✓ Random generation works")
    print(f"  Generated {nbc.num_codewords} codewords")


def test_library_indices():
    """Test library index packing mode."""
    print("\nTest 12: Library index packing")
    
    from dnabyte.library import Library
    import os
    
    # Skip if library file doesn't exist
    lib_path = './tests/testlibraries/20bp_Lib.csv'
    if not os.path.exists(lib_path):
        print("⊘ Library index test skipped (20bp_Lib.csv not found)")
        return
    
    # Load library
    lib = Library(structure='linear_assembly', filename=lib_path)
    
    # Create data with indices
    data = [
        [[0, 1, 2], [3, 4, 5]],
        [[6, 7, 8]]
    ]
    
    nbc = NucleobaseCode.from_library_indices(data, lib)
    
    assert nbc.pack_mode == 'library_indices'
    assert nbc.num_codewords == 2
    
    # Unpack should give library sequences
    unpacked = nbc.unpack()
    assert isinstance(unpacked[0][0][0], str)
    assert len(unpacked[0][0]) == 3
    
    print("✓ Library index packing works")
    print(f"  Index bits: {nbc.index_bits}")
    print(f"  Packed size: {len(nbc._packed_data)} bits")
    print(f"  _structure: {nbc._structure}")
    print(f"  _packed_data: {nbc._packed_data.to01()}")


def test_string_representation():
    """Test __str__ method."""
    print("\nTest 13: String representation")
    
    data = [['ATCG', 'GGCC']]
    nbc = NucleobaseCode(data)
    
    str_repr = str(nbc)
    assert "NucleobaseCode" in str_repr
    assert "packed=True" in str_repr
    assert "mode=dna" in str_repr
    
    print("✓ String representation works")


def run_all_tests():
    """Run all test functions."""
    print("Running NucleobaseCode v2 tests...")
    print("=" * 70)
    
    test_functions = [
        test_basic_init,
        test_library_indices_init,
        test_single_bitarray_storage,
        test_nested_structure_preservation,
        test_2bit_encoding,
        test_indexing,
        test_iteration,
        test_validation,
        test_already_packed,
        test_memory_efficiency,
        test_complement,
        test_random_generation,
        test_string_representation,
    ]
    
    passed = 0
    failed = 0
    
    for test_func in test_functions:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"✗ {test_func.__name__} failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("=" * 70)
    print(f"Tests passed: {passed}/{len(test_functions)}")
    
    if failed > 0:
        print(f"Tests failed: {failed}")
        sys.exit(1)
    else:
        print("All tests passed!")
        sys.exit(0)


if __name__ == '__main__':
    run_all_tests()
