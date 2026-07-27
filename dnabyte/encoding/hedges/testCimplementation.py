from hedges_python import *

original = b"Hello HEDGES"

dna = encode_hedges(original)

print(len(dna))
print(dna[0])

decoded = decode_hedges(dna)

print(bytes(decoded))