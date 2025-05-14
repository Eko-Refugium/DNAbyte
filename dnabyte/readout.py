import random
import re

def ErrorfunctionLenght(length):
    return min((length-1)/10,1)
    #return 0.1

def ErrorfunctionMagnitude(length):
    return random.poisson(lam=1)
    #return 1

# def DNA_to_quaternary(DNA_string):
#     mapping = {'A': '0', 'T': '1', 'C': '2', 'G': '3'}
#     return ''.join(mapping[letter] for letter in DNA_string)

def DNA_to_quaternary(DNA_string):
    mapping = {'A': '0', 'T': '1', 'C': '2', 'G': '3'}
    if isinstance(DNA_string, list):
        DNA_string = DNA_string[0]
    return ''.join(mapping[base] for base in DNA_string)

def quaternary_to_DNA(quaternarystring):
    mapping = {'0': 'A', '1': 'T', '2': 'C', '3': 'G'}
    return ''.join(mapping[digit] for digit in quaternarystring)

def find_and_count_sequential_digits(s):
    matches = re.findall(r'(\d)\1*', s)
    return [(m.start(), len(m.group())) for m in re.finditer(r'(\d)\1*', s)]


def replace_random_base4(s, start, length, FunctionOfErrorDependingOnLength, FunctionOfMagnitudeOfError, q):
    base4_numbers = range(q)
    substring = s[start:start+length]
    
    # Probability of replacement is proportional to the length of the substring
    if random.random() < FunctionOfErrorDependingOnLength(length):
        # Number of replacements is also proportional to the length of the substring
        num_replacements = random.randint(0,int(FunctionOfMagnitudeOfError(length)))
        for _ in range(num_replacements):
            index_to_replace = random.randint(0, len(substring)-1)
            new_number = random.choice(base4_numbers)
            substring = substring[:index_to_replace] + str(new_number) + substring[index_to_replace+1:]
    return s[:start] + substring + s[start+length:]

def Implement_homopolymer_Error(DNA_string,base,FunctionOfErrorDependingOnLength, FunctionOfMagnitudeOfError):
    NumberString=DNA_to_quaternary(DNA_string)
    WhereandHowlong=find_and_count_sequential_digits(NumberString)
    for i in range(len(WhereandHowlong)):
        NumberString=replace_random_base4(NumberString, WhereandHowlong[i][0], WhereandHowlong[i][1], FunctionOfErrorDependingOnLength, FunctionOfMagnitudeOfError, base)
    
    ErrDNA_string=quaternary_to_DNA(NumberString)
    return ErrDNA_string


def sequencing_simulator(stored_data):
    sequenced_data = Implement_homopolymer_Error(stored_data, 3, ErrorfunctionLenght, ErrorfunctionMagnitude)
    return sequenced_data


