import dnabyte.error_correction.fountain_code.fountain_code_encoding_nit as fce
import dnabyte.error_correction.fountain_code.fountain_code_decoding_nit as fcd
import numpy as np
import random

class Symbol:
    def __init__(self, index, degree, data, neighbors):
        self.index = index
        self.degree = degree
        self.data = data
        self.neighbors=neighbors

# Function to find a random
# number between 0 or 1
def findRandom():
     
    # Generate the random number
    num = random.randint(0, 1)
 
    # Return the generated number
    return num
 
# Function to generate a random
# binary string of length N
def generateBinaryString(N):
     
    # Stores the empty string
    S = ""
 
    # Iterate over the range [0, N - 1]
    for i in range(N):
         
        # Store the random number
        x = findRandom()
 
        # Append it to the string
        S += str(x)
     
    return S

def translate_dna(dna_string):
    # Create a mapping from DNA bases to numbers
    dna_mapping = {'A': '0', 'C': '1', 'G': '0', 'T': '1'}

    # Translate the DNA string into numbers
    binary_string = ''.join(dna_mapping[dna] for dna in dna_string)

    # Represent the binary string as a binary number
    binary_number = int(binary_string, 2)

    
    base10 = np.base_repr(binary_number, base=10)
   

    return base10

def translate_binary(binary_string):
    # Create a mapping from binary digits to DNA bases
    binary_mapping = {0: ['A', 'G'], 1: ['C', 'T']}

    # Translate the binary string into DNA bases
    dna_string = ''.join(binary_mapping[int(digit)][i % 2] for i, digit in enumerate(binary_string))

    return dna_string

def initilation(sortsym):
    finallist=[]
    while sortsym[0].degree == 1:
        finallist.append(sortsym[0])
        sortsym.pop(0)
    return finallist,sortsym

def reprodice(finalist,restlist,finalnumber):
    check=0
    while len(finalist)<finalnumber:
        check=check+1
        
        for j in finalist:
            i=0
            upperlimit=len(restlist)
            while i < upperlimit:
            
                if j.neighbors[0] in restlist[i].neighbors:
                    check=0
                    restlist[i].data = np.bitwise_xor(restlist[i].data, j.data)
                    restlist[i].degree = restlist[i].degree - 1 
                    restlist[i].neighbors.remove(j.neighbors[0])
                    finaldata=[]
                    for k in finalist:
                        finaldata.append(k.data)

                    if restlist[i].degree == 1 and restlist[i].data not in finaldata:
                        finalist.append(restlist[i])
                        del restlist[i]
                        upperlimit=upperlimit-1
                        i=i-1
                    if restlist[i].degree == 1 and restlist[i].data in finaldata:
                        
                        del restlist[i]
                        upperlimit=upperlimit-1
                        i=i-1
                i=i+1
        
        if check>100:
            break
    
    return finalist

def translate_binary(binary_string):
    # Create a mapping from binary digits to DNA bases
    binary_mapping = {0: ['A', 'G'], 1: ['C', 'T']}

    # Translate the binary string into DNA bases
    dna_string = ''.join(binary_mapping[int(digit)][i % 2] for i, digit in enumerate(binary_string))

    return dna_string
    

# Convert the numbers to binary
binary_numbers = []
for i in range(100):
    binary_numbers.append(generateBinaryString(10))


dna_sequences = []
for i in binary_numbers:
    
    dna_sequences.append(translate_binary(i))
dna_sequences=sorted(list(set(dna_sequences)))

test=fce.dnafountaincode(dna_sequences, 150)
for i in test:
    i.data=int(translate_dna(i.data))

decoded=fcd.fountaindecode(test,100,10)


