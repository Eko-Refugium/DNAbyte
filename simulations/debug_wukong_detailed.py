#!/usr/bin/env python3
"""Debug Wukong encoding/decoding at 0.01 error rate with detailed tracing."""

from dnabyte.params import Params
from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.data_classes import Data
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus

def default_wukong():
    return {
        'name': 'wukong_default',
        'sequence_length': 200,
        'encoding_method': 'wukong',
        'max_homopolymer': 3,
        'min_gc_content': 0.4,
        'max_gc_content': 0.6,
        'rule_num': 1,
        'rs_num': 0,
        'add_redundancy': False,
        'add_primer': False,
    }

# Create initial params
base_params = {
    'filename': 'testfilesimsall.txt',
    'binarization_method': 'default',
    'synthesis_method': 'nosynthpoly',
    'sequencing_method': 'iid',
    'recovery_method': 'simple',
    'min_coverage': 1,
    'clustering_method': 'kmere_cluster',
    'storage_conditions': None,
    'kmer_seed': 42,
    'mean': 1,
    'std_dev': 0,
    "iid_error_rate": 0.0,
    "iid_substitution_rate": 0.0,
    "iid_insertion_rate": 0.0,
    "iid_deletion_rate": 0.0,
}

wukong_params = base_params.copy()
wukong_params.update(default_wukong())

print("=" * 80)
print("WUKONG DEBUG: Testing with 0.0 error rate first")
print("=" * 80)

params = Params(**wukong_params)

file_paths = ['./tests/testfiles/' + params.filename]
data_obj = Data(file_paths=file_paths)

binarizer = Binarize(params)
binary_code = binarizer.binarize(data_obj)
print(f"Original binary code length: {len(binary_code.data[0])} bits")

coder = Encode(params)
data_enc, info = coder.encode(binary_code)
print(f"Encoded sequences: {len(data_enc.data)} sequences")
print(f"First 5 sequences: {[seq[:50] for seq in data_enc.data[:5]]}")

syn = SimulateSynthesis(params)
data_syn, info = syn.simulate(data_enc)
print(f"After synthesis: {len(data_syn.data)} sequences (mean={params.mean})")

sto = SimulateStorage(params)
data_sto, info = sto.simulate(data_syn)
print(f"After storage: {len(data_sto.data)} sequences")

seq = SimulateSequencing(params)
data_seq, info = seq.simulate(data_sto)
print(f"After sequencing (0.0 error): {len(data_seq.data)} sequences")

# Try clustering + recovery
cluster_obj = Cluster(params)
data_cluster, info = cluster_obj.cluster(data_seq)
print(f"After clustering: {len(data_cluster.data)} clusters")

consensus_obj = Consensus(params)
data_cor, info = consensus_obj.call(data_cluster)
print(f"After consensus: {len(data_cor.data)} sequences")

data_dec, valid, info = coder.decode(data_cor)
print(f"Decoded: valid={valid}, decoded length={len(data_dec.data[0]) if data_dec.data else 0} bits")

if valid:
    comparison, res = data_dec.compare(data_dec, binary_code)
    print(f"[OK] 0.0 error rate: {comparison}")
else:
    print(f"[FAIL] 0.0 error rate: Decode failed")

print("\n" + "=" * 80)
print("WUKONG DEBUG: Testing with 0.01 error rate")
print("=" * 80)

# Now test with 0.01 error rate
wukong_params['iid_error_rate'] = 0.01
wukong_params['iid_substitution_rate'] = 0.01
wukong_params['iid_insertion_rate'] = 0.0
wukong_params['iid_deletion_rate'] = 0.0

params = Params(**wukong_params)

file_paths = ['./tests/testfiles/' + params.filename]
data_obj = Data(file_paths=file_paths)

binarizer = Binarize(params)
binary_code = binarizer.binarize(data_obj)
print(f"Original binary code length: {len(binary_code.data[0])} bits")

coder = Encode(params)
data_enc, info = coder.encode(binary_code)
print(f"Encoded sequences: {len(data_enc.data)} sequences")
print(f"First sequence length: {len(data_enc.data[0])}")

syn = SimulateSynthesis(params)
data_syn, info = syn.simulate(data_enc)
print(f"After synthesis: {len(data_syn.data)} sequences")

sto = SimulateStorage(params)
data_sto, info = sto.simulate(data_syn)
print(f"After storage: {len(data_sto.data)} sequences")

seq = SimulateSequencing(params)
data_seq, info = seq.simulate(data_sto)
print(f"After sequencing (0.01 error): {len(data_seq.data)} sequences")
print(f"First 3 sequences after errors:")
for i, s in enumerate(data_seq.data[:3]):
    print(f"  [{i}] len={len(s)}: {s[:60]}...")

# Try clustering + recovery
print(f"\nClustering with kmere_cluster (k={params.kmer_size_cluster})...")
cluster_obj = Cluster(params)
data_cluster, info = cluster_obj.cluster(data_seq)
print(f"After clustering: {len(data_cluster.data)} clusters")
print(f"Cluster sizes: {[len(cluster) for cluster in data_cluster.data[:5]]}")

consensus_obj = Consensus(params)
data_cor, info = consensus_obj.call(data_cluster)
print(f"After consensus: {len(data_cor.data)} sequences")
print(f"First consensus sequence length: {len(data_cor.data[0]) if data_cor.data else 'N/A'}")

data_dec, valid, info = coder.decode(data_cor)
print(f"Decoded: valid={valid}, decoded length={len(data_dec.data[0]) if data_dec.data else 0} bits")

if valid:
    comparison, res = data_dec.compare(data_dec, binary_code)
    print(f"[OK] 0.01 error rate: {comparison}")
else:
    print(f"[FAIL] 0.01 error rate: Decode failed")
    if data_dec.data and len(data_dec.data) > 0:
        print(f"   Decoded {len(data_dec.data[0])} bits (expected {len(binary_code.data[0])})")
