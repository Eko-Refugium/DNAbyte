import math as m
codeword_length = 500

def decimal_to_binary(n): 
    return bin(n).replace("0b", "") 

thebitsforzfill = m.ceil(m.log2(codeword_length))

print('thebitesforzfill:', thebitsforzfill)

codeword = '00101101011111010100001010001010001'

print('codeword:', codeword)

codeword = str(decimal_to_binary(0)).zfill(thebitsforzfill * 2) + codeword

print('codeword_after:', codeword)