import os
from hedges_python import *


cfg = HedgesConfig()


for size in [
    1,
    5,
    13,
    861,
    875,
    2500,
    6885,
    7000,
    20000,
    100000
]:

    # create deterministic binary pattern
    original = ''.join(
        ['1' if (i % 2 == 0) else '0' for i in range(size)]
    )


    # convert bits -> bytes
    byte_data = int(original, 2).to_bytes(
        (len(original)+7)//8,
        byteorder="big"
    )


    # encode
    dna = encode_hedges(
        byte_data,
        cfg
    )


    # decode
    decoded = decode_hedges(
        dna,
        cfg
    )


    # bytes -> bits
    decoded_binary = ''.join(
        f'{b:08b}' for b in decoded
    )


    # only compare the original amount of bits
    decoded_bits = decoded_binary[-size:]


    wherewrong = [
        i for i in range(
            min(size, len(decoded_bits))
        )
        if original[i] != decoded_bits[i]
    ]


    print(
        size,
        len(dna),
        decoded_bits == original,
        "decoded bits:",
        len(decoded_bits),
        "wrong:",
        len(wherewrong)
    )

    if wherewrong:
        print("first error:", wherewrong[0])

import os
from hedges_python import *


cfg = HedgesConfig()

# change HEDGES parameters
cfg.coderate = 4

cfg.total_strand_length = 250

cfg.strand_id_bytes = 2
cfg.strand_runout_bytes = 2

cfg.gc_window = 14
cfg.gc_max = 9

cfg.max_homopolymer = 3

cfg.heap_limit = 500000


# primers
cfg.left_primer = "TCGAAGTCAGCGTGTATTGTATG"
cfg.right_primer = "TAGTGAGTGCGATTAAGCGTGTT"

for size in [
    1,
    5,
    13,
    861,
    875,
    2500,
    6885,
    7000,
    20000,
    100000
]:

    # create deterministic binary pattern
    original = ''.join(
        ['1' if (i % 2 == 0) else '0' for i in range(size)]
    )


    # convert bits -> bytes
    byte_data = int(original, 2).to_bytes(
        (len(original)+7)//8,
        byteorder="big"
    )


    # encode
    dna = encode_hedges(
        byte_data,
        cfg
    )


    # decode
    decoded = decode_hedges(
        dna,
        cfg
    )


    # bytes -> bits
    decoded_binary = ''.join(
        f'{b:08b}' for b in decoded
    )


    # only compare the original amount of bits
    decoded_bits = decoded_binary[-size:]


    wherewrong = [
        i for i in range(
            min(size, len(decoded_bits))
        )
        if original[i] != decoded_bits[i]
    ]


    print(
        size,
        len(dna),
        decoded_bits == original,
        "decoded bits:",
        len(decoded_bits),
        "wrong:",
        len(wherewrong)
    )

    if wherewrong:
        print("first error:", wherewrong[0])

# import os
# from hedges_python import *

# for size in [
#     1,
#     40,
#     100,
#     6885,
#     7000,
#     20000,
#     100000
# ]:
#     original = ''.join(['1' if (i % 2 == 0) else '0' for i in range(size)])

#     byte_data = int(original, 2).to_bytes(
#                 (len(original)+7)//8,
#                 byteorder="big"
#             )

#     dna = encode_hedges(byte_data)
#     decoded = decode_hedges(dna)
    
#     decoded_binary = ''.join(f'{b:08b}' for b in decoded)

#     wherewrong = [i for i in range(min(size, len(decoded_binary[-size:]))) if original[i] != decoded_binary[-size:][i]]

#     print(
#         size,
#         len(dna),
#         decoded_binary[-size:] == original,
#         len(decoded_binary[-size:]),
#         len(original),
#         # print(wherewrong)
#     )

# original = "00110010001111110100001100111111011101000111011001000010001010100100000101011110010110110101000001110110001001010101011100111011001100110110100101100011010101000011010101100111001100000111000001000100001000110011110101011010001111110110001001100011011101110011001101001000011010010010000100100110010000100100110101011001"

# byte_data = int(original, 2).to_bytes(
#             (len(original)+7)//8,
#             byteorder="big"
#         )

# dna = encode_hedges(byte_data)

# print("Number of strands:", len(dna))
# print("First strand length:", len(dna[0]))
# print("First strand:", dna[0])
# print("Second strand:", dna[1])

# decoded = decode_hedges(dna)

# decoded_binary = ''.join(f'{b:08b}' for b in decoded)

# # print(decoded_binary)
# print(decoded_binary == original)#