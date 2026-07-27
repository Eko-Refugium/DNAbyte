from hedges_python import *

original = "00110010001111110100001100111111011101000111011001000010001010100100000101011110010110110101000001110110001001010101011100111011001100110110100101100011010101000011010101100111001100000111000001000100001000110011110101011010001111110110001001100011011101110011001101001000011010010010000100100110010000100100110101011001"

byte_data = int(original, 2).to_bytes(
            (len(original)+7)//8,
            byteorder="big"
        )

dna = encode_hedges(byte_data)

print("Number of strands:", len(dna))
print("First strand length:", len(dna[0]))
print("First strand:", dna[0])
print("Second strand:", dna[1])

decoded = decode_hedges(dna)

decoded_binary = ''.join(f'{b:08b}' for b in decoded)

# print(decoded_binary)
print(decoded_binary == original)