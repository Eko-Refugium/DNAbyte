#!/usr/bin/env python
"""Debug script to test clustering + recovery with duplicates"""

from dnabyte.post_sequencing.clustering.primer_grouper.clustering import PrimerGrouper
from dnabyte.post_sequencing.recovery.debruijn.recovery import ConsensusBuilder

# Use actual church-encoded sequences, duplicated 10 times (like synthesis output)
seq1 = 'GGAAACACGCCTTCACCCAGAACAGGAACGAAGGAACGAAGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACAACGGAACGAAGGAACAACGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACCGAAGGCACCCTGGTACTT'
seq2 = 'GGAAACACGCCTTCACCCAGAGAACAGGAACAACGGAACAACGGAACAACGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACAACGGAACGAAGGAACGAAGGAACGAAGGCCGAAGGCACCCTGGTACTT'

seqs_with_coverage = [seq1] * 10 + [seq2] * 10

print("="*80)
print("CLUSTERING WITH DUPLICATES")
print("="*80)

grouper = PrimerGrouper()

print(f"\nInput: {len(seqs_with_coverage)} sequences")
print(f"  10 copies of seq1 + 10 copies of seq2")
print(f"  All 200bp each")

grouped = grouper.cluster(seqs_with_coverage)

print(f"\nGrouped into {len(grouped)} groups:")
for key, payloads in grouped.items():
    left, right, *rest = key
    hash_val = rest[0] if rest else "N/A"
    print(f"  Key: ({left[:10]}..., ...{right[-10:]}, {hash_val})")
    print(f"    Payloads: {len(payloads)} items, each {len(payloads[0])}bp")
    if payloads[0] != payloads[1]:
        print(f"    WARNING: Payloads are not identical!")
        print(f"      [0]: {payloads[0][:40]}...")
        print(f"      [1]: {payloads[1][:40]}...")

print("\n" + "="*80)
print("RECOVERY")
print("="*80)

recovery = ConsensusBuilder(method='debruijn', kmer_size=21, min_coverage=1)

print(f"\nInput from clustering: {len(grouped)} groups")
result = recovery.recover(grouped)
consensus_seqs, stats = result

print(f"\nResult: {len(consensus_seqs)} consensus sequences")
for i, seq in enumerate(consensus_seqs):
    print(f"  [{i}] Length: {len(seq)}bp (expected 200bp)")
    print(f"      {seq[:50]}...{seq[-30:]}")

print(f"\nStats: {stats}")
