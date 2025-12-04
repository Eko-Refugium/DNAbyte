import csv
import random
from pyxdameraulevenshtein import damerau_levenshtein_distance
import itertools
import matplotlib.pyplot as plt
import numpy as np
import math as m


def Complimentmap(string):
    compliment = ''
    for i in string:
        if i == 'A':
            compliment += 'T'
        if i == 'T':
            compliment += 'A'
        if i == 'C':
            compliment += 'G'
        if i == 'G':
            compliment += 'C'
    compliment = compliment[::-1]
    return compliment


def Checkingcompliment(element1,element2,n):
    le1=element1[:n]
    re1=element1[n:]
    le2=element2[:n]
    re2=element2[n:]
    if Complimentmap(le1) == le2 or Complimentmap(re1) == re2:
        return True
    else:
        return False
    
def generate_balanced_string(size):
    while True:
        string = ''.join(random.choice('ACTG') for _ in range(size))
        at_count = string.count('A') + string.count('T')
        cg_count = string.count('C') + string.count('G')
        if abs(at_count - cg_count) <= 0:  # Allow a difference of up to 2 for "somewhat balanced"
            return string
        
def generate_balanced_string_nohomo(size):
    while True:
        string = ''.join(random.choice('ACTG') for _ in range(size))
        at_count = string.count('A') + string.count('T')
        cg_count = string.count('C') + string.count('G')
        
        # Check for balanced AT and CG counts
        if abs(at_count - cg_count) <= 3:
            # Check for consecutive identical characters
            if all(string[i] != string[i + 1] for i in range(len(string) - 1)):
                return string

def similaritycheck(string1,string2):
    count=0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            count+=1
    return count

def generate_random_string(n):
    return ''.join(random.choice('ATCG') for _ in range(n))

def combine_strings(list1, list2):
    return [s1 + s2 for s1 in list1 for s2 in list2]

def moving_window(sequence, window_size):
    
    for i in range(len(sequence) - window_size + 1):
        yield sequence[i:i + window_size]

def generate_substrings(s):
    substrings = []
    
    # Generate all prefixes
    for i in range(1, len(s) + 1):
        substrings.append(s[:i])
    
    # Generate all suffixes
    for i in range(1, len(s)):
        substrings.append(s[i:])
    
    return substrings


def generate_positional_library(lengthmessage, amountofmessages, lengthgeneric, lengthposition, amountofposition):
    
    amoauntofgeneric = 2
    genericset = []
    messageset = []
    positionset = []


    tester=0
    while len(genericset) < amoauntofgeneric:
        if tester == 0:
            genericset.append(generate_balanced_string_nohomo(lengthgeneric))
            
            
            tester+=1
            
        newgeneric = generate_balanced_string_nohomo(lengthgeneric)
        if damerau_levenshtein_distance(genericset[0],newgeneric) < lengthgeneric/2:
            genericset.append(newgeneric)



    genericorder =[]
    reversestrands = []
    tester=0

    while len(messageset) <= amountofmessages:
        # print(len(messageset)<amountofmessages,'ifcon')
        if tester == 0:
            messageset.append(generate_balanced_string(lengthmessage))
            tester+=1

        newmessage = generate_random_string(lengthmessage)
        for alredyin in messageset:
            if damerau_levenshtein_distance(alredyin,newmessage) < lengthmessage/2:
                break
        else:
            genericordertemp=genericset[0]+newmessage+genericset[1]
            messagereversestrand = Complimentmap(newmessage)
            counter=0
            for window in moving_window(genericordertemp, lengthmessage):
                counter+=1
                if counter not in range(lengthgeneric-5,lengthgeneric+lengthmessage+5) and damerau_levenshtein_distance(window,newmessage) < 5:
                    
                    break
            else:
                genericorder.append(genericordertemp)
                reversestrands.append(messagereversestrand)
                messageset.append(newmessage)
        # print(len(messageset),amountofmessages,'ifcon')
   

    # tester=0

    while len(positionset) < amountofposition*2:
        # if tester == 0:
        #     positionset.append(generate_balanced_string(lengthposition))
        #     tester+=1

        newposition = generate_random_string(lengthposition)
        for alredyin in positionset:
            if damerau_levenshtein_distance(alredyin,newposition) < lengthposition/2:
                break
        else:
            positionset.append(Complimentmap(newposition))

    file=open(f'./tests/testlibraries/lib_positional_{lengthmessage}({amountofmessages})m_{lengthgeneric}({amoauntofgeneric})g_{lengthposition}({amountofposition})p.csv', 'w')
    file.write('Messages\n')
    for i in range(len(messageset)):
        file.write(messageset[i]+'\n')
    file.write('Generic\n')
    for i in range(len(genericset)):
        file.write(genericset[i]+'\n')
    file.write('Connector\n')
    for i in range(len(positionset)):
        file.write(positionset[i]+'\n')
    file.close()

            

