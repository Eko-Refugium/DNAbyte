#!/usr/bin/env python
"""Debug script to trace De Bruijn graph construction"""

from dnabyte.post_sequencing.recovery.debruijn.recovery import (
    build_debruijn_graph,
    find_consensus_path_greedy,
    extract_consensus_from_sequences
)

# Single 160bp payload (church encrypted middle part)
payload = 'AACAGGAACGAAGGAACGAAGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACAACGGAACGAAGGAACAACGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAAC'

print("="*80)
print("PAYLOAD ANALYSIS")
print("="*80)
print(f"\nPayload length: {len(payload)}bp")
print(f"Content: {payload[:50]}...{payload[-50:]}\n")

# Test with different kmer sizes
for k in [11, 15, 21, 25]:
    print(f"\n{'='*80}")
    print(f"KMER SIZE: {k}")
    print(f"{'='*80}")
    
    # Test 1: Single sequence
    print(f"\n1. Single sequence")
    try:
        result = extract_consensus_from_sequences([payload], kmer_size=k, min_coverage=1)
        print(f"   Result length: {len(result)}bp (expected {len(payload)}bp)")
        print(f"   Content: {result[:50]}...{result[-30:]}")
        if len(result) != len(payload):
            print(f"   ERROR: Length mismatch!")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Test 2: 10 identical sequences (coverage)
    print(f"\n2. 10 identical sequences (coverage)")
    try:
        result = extract_consensus_from_sequences([payload] * 10, kmer_size=k, min_coverage=1)
        print(f"   Result length: {len(result)}bp (expected {len(payload)}bp)")
        print(f"   Content: {result[:50]}...{result[-30:]}")
        if len(result) != len(payload):
            print(f"   ERROR: Length mismatch! Got {len(result)}, expected {len(payload)}")
    except Exception as e:
        print(f"   Error: {e}")

# Test graph structure
print(f"\n{'='*80}")
print(f"GRAPH STRUCTURE ANALYSIS (k=21)")
print(f"{'='*80}")

k = 21
graph = build_debruijn_graph([payload], kmer_size=k, min_coverage=1)
print(f"\nGraph stats:")
print(f"  Total nodes: {len(graph.graph)}")
print(f"  Total kmers: {sum(len(edges) for edges in graph.graph.values())}")
print(f"  Start nodes: {graph.get_start_nodes()}")
print(f"  End nodes: {graph.get_end_nodes()}")

# Try to find path
path = find_consensus_path_greedy(graph)
print(f"\nGreedy traversal result:")
print(f"  Path length: {len(path)}bp (expected {len(payload)}bp)")
if len(path) != len(payload):
    print(f"  ERROR: Path length mismatch!")
