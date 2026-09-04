#!/usr/bin/env python3
"""Minimal reproducer of Cost_analysis pipeline for wukong at 0.01 error"""

from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.data_classes import Data
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus

# Wukong parameters after error correction tuning (from encode_and_count_with_error_correction)
params_dict = {
    'filename': 'testfilesimsall.txt',
    'binarization_method': 'default',
    'encoding_method': 'wukong',
    'sequence_length': 200,
    'max_homopolymer': 3,
    'min_gc_content': 0.4,
    'max_gc_content': 0.6,
    'rule_num': 1,
    'rs_num': 2,  # After tuning
    'add_redundancy': True,  # After tuning
    'add_primer': False,
    'synthesis_method': 'nosynthpoly',
    'sequencing_method': 'iid',
    'recovery_method': 'simple',
    'min_coverage': 1,
    'clustering_method': 'kmere_cluster',
    'storage_conditions': None,
    'kmer_seed': 42,
    'mean': 3,  # After tuning
    'std_dev': 0,
    "iid_error_rate": 0.01,
    "iid_substitution_rate": 0.01,
    "iid_insertion_rate": 0.0,
    "iid_deletion_rate": 0.0,
}

print("=" * 80)
print("Minimal Wukong Pipeline at 0.01 Error Rate")
print("=" * 80)

params = Params(**params_dict)

# Step 1: Load and binarize
print("\n1. Binarizing...")
file_paths = ['./tests/testfiles/' + params.filename]
data_obj = Data(file_paths=file_paths)
binarizer = Binarize(params)
binary_code = binarizer.binarize(data_obj)

print(f"   Binary code length: {len(binary_code.data)} bits")
print(f"   First 100 chars: {binary_code.data[:100]}")

# Step 2: Encode
print("\n2. Encoding with Wukong (rs_num=2, add_redundancy=True)...")
coder = Encode(params)
data_enc, info = coder.encode(binary_code)
print(f"   Encoded sequences: {len(data_enc.data)}")
print(f"   First sequence length: {len(data_enc.data[0]) if data_enc.data else 'N/A'}")

# Step 3: Synthesize (mean=3 copies)
print("\n3. Synthesizing (mean=3)...")
syn = SimulateSynthesis(params)
data_syn, info = syn.simulate(data_enc)
print(f"   After synthesis: {len(data_syn.data)} sequences (expected ~{len(data_enc.data)*3})")

# Step 4: Storage
print("\n4. Storing...")
sto = SimulateStorage(params)
data_sto, info = sto.simulate(data_syn)
print(f"   After storage: {len(data_sto.data)} sequences")

# Step 5: Sequencing with 0.01 error rate
print("\n5. Sequencing with 0.01 error rate...")
seq = SimulateSequencing(params)
data_seq, info = seq.simulate(data_sto)
print(f"   After sequencing: {len(data_seq.data)} sequences")

# Step 6: Clustering
print("\n6. Clustering (kmere_cluster)...")
cluster_obj = Cluster(params)
data_cluster, info = cluster_obj.cluster(data_seq)
print(f"   After clustering: {len(data_cluster.data)} clusters/groups")
if hasattr(data_cluster, 'data') and isinstance(data_cluster.data, dict):
    print(f"   Cluster sizes: {[len(v) for v in list(data_cluster.data.values())[:5]]}")

# Step 7: Consensus/Recovery
print("\n7. Consensus (simple majority voting)...")
consensus_obj = Consensus(params)
data_cor, info = consensus_obj.call(data_cluster)
print(f"   After consensus: {len(data_cor.data)} sequences")
if data_cor.data:
    print(f"   First sequence length: {len(data_cor.data[0]) if isinstance(data_cor.data, list) else 'N/A'}")

# Step 8: Decode
print("\n8. Decoding...")
data_dec, valid, info = coder.decode(data_cor)
print(f"   Decoded: valid={valid}, decoded length={len(data_dec.data)} bits")

# Step 9: Compare
print("\n9. Comparing...")
comparison, res = data_dec.compare(data_dec, binary_code)
print(f"   Result: {comparison}")
if comparison == 'SUCCESS':
    print("   [OK] 0.01 error rate recovery successful!")
else:
    print("   [FAIL] 0.01 error rate recovery failed")
    print(f"   Original: {len(binary_code.data)} bits")
    print(f"   Decoded:  {len(data_dec.data)} bits")
