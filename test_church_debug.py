#!/usr/bin/env python
"""Debug script to trace church encoding pipeline"""

from dnabyte.params import Params
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.encode import Encode
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.process import Process

# Setup parameters exactly like test
params = Params(
    name='end2end_church_basic',
    filename='textfile_40b.txt',
    
    # encoding parameters
    encoding_method='church',
    binarization_method='default',
    sequence_length=200,
    max_homopolymer=2,
    rs_num=0,
    add_redundancy=True,
    add_primer=True,
    primer_length=20,
    
    # error channels
    storage_conditions=None,
    synthesis_method='nosynthpoly',
    clustering_method='primer_grouper',
    recovery_method='debruijn'
)

# Set defaults for synthesis
if not hasattr(params, 'mean'):
    params.mean = 1
if not hasattr(params, 'std_dev'):
    params.std_dev = 0

# Create simple test binary
binary_code = BinaryCode("00110001001100010011000100110000001100010011000000110001001100000011000100110000001100000011000100110000001100000011000100110000001100010011000000110001001100000011000000110000001100000011000000110001001100000011000100110000001100010011000000110001001100000011000100110000001100010011000000110001001100010011000100110001001100010011000000110001001100000011000100110000001100010011000000110000001100000011000100110000001100010011000000110001")

# STEP 1: Encode (use wrapper)
print("\n=== STEP 1: ENCODE ===")
from dnabyte.encode import Encode
enc = Encode(params)
data_enc, info = enc.encode(binary_code)
print(f"Result type: {type(data_enc)}")
print(f"Result.data type: {type(data_enc.data)}")
print(f"Number of sequences: {len(data_enc.data)}")
print(f"First sequence length: {len(data_enc.data[0])}")
print(f"Params has metadata: {hasattr(params, 'church_fasta_metadata')}")
print(f"Metadata: {params.church_fasta_metadata[:80]}")

# STEP 2: Synthesis
print("\n=== STEP 2: SYNTHESIS ===")
syn = SimulateSynthesis(params)
data_syn, info = syn.simulate(data_enc)
print(f"Result type: {type(data_syn)}")
print(f"Result.data type: {type(data_syn.data)}")
print(f"Number of sequences: {len(data_syn.data)}")
print(f"First sequence length: {len(data_syn.data[0])}")
print(f"Params has metadata: {hasattr(params, 'church_fasta_metadata')}")

# STEP 3: Convert to InSilicoDNA (if needed)
print("\n=== STEP 3: ENSURE INSILICODNA ===")
if not isinstance(data_syn, InSilicoDNA):
    print("Converting to InSilicoDNA...")
    data_syn = InSilicoDNA(data_syn.data)
else:
    print("Already InSilicoDNA")
print(f"Result type: {type(data_syn)}")
print(f"Number of sequences: {len(data_syn.data)}")
print(f"First sequence length: {len(data_syn.data[0])}")

# STEP 4: Process (clustering+recovery)
print("\n=== STEP 4: PROCESS ===")
processor = Process(params)
data_cor, info = processor.process(data_syn)
print(f"Result type: {type(data_cor)}")
print(f"Result.data type: {type(data_cor.data)}")
print(f"Number of sequences: {len(data_cor.data)}")
for i, seq in enumerate(data_cor.data[:2]):
    print(f"Seq {i}: length={len(seq)}, first 50bp={seq[:50]}, last 30bp={seq[-30:]}")
print(f"Params has metadata: {hasattr(params, 'church_fasta_metadata')}")

# STEP 5: Decode
print("\n=== STEP 5: DECODE ===")
print(f"Data to decode: {len(data_cor.data)} sequences")
print(f"First sequence: {data_cor.data[0][:80]}...")
decoded, valid, info = enc.decode(data_cor)
print(f"Decoded type: {type(decoded)}")
print(f"Decoded value: {str(decoded.data)[:100] if hasattr(decoded, 'data') else str(decoded)[:100]}")
print(f"Decoded length: {len(decoded.data) if hasattr(decoded, 'data') else len(decoded) if isinstance(decoded, str) else 'N/A'}")
print(f"Valid: {valid}")
