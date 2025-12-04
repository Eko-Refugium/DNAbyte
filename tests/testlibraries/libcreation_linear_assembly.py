
import csv
import random
import itertools
import matplotlib.pyplot as plt
import numpy as np
import math


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
    
def generate_combinations(n, current='', combinations=[]):
    if n == 0:
        if (current[0] in 'AC') and (current[-1] in 'TG'):
            combinations.append(current)
    else:
        for letter in 'ATCG':
            if not current or current[-1] != letter:
                generate_combinations(n-1, current+letter, combinations)
    return combinations

def similaritycheck(string1,string2):
    count=0
    for i in range(len(string1)):
        if string1[i] != string2[i]:
            count+=1
    return count

def generate_balanced_string(size,absdiff):
    while True:
        string = ''.join(random.choice('ACTG') for _ in range(size))
        at_count = string.count('A') + string.count('T')
        cg_count = string.count('C') + string.count('G')
        if abs(at_count - cg_count) <= absdiff:  # Allow a difference of up to 2 for "somewhat balanced"
            return string
def simplecration(oligosize,motiveamounts):
    
    

    #data=[]
    # with open('24oligolibrand.txt') as csv_file:
    #     csv_reader = csv.reader(csv_file, delimiter=',')
    #     line_count = 0
    #     for row in csv_reader:
    #         data.append(row[0])
    #         line_count+=1
    # #data = generate_combinations(10)
    amountofoligos=(2*motiveamounts)**2

    counterenergy=0

    runthroughsize=oligosize
    runthroughsizehalf=int(runthroughsize/2)
    erngylimit=-counterenergy*3.5/2-5.3
    similaritycounter=math.floor(runthroughsize/6-1)
    absdiff=math.floor(runthroughsize/8)

    leftlist=[]
    countleft=0
    while countleft<motiveamounts:
        selected=generate_balanced_string(runthroughsizehalf,absdiff)
        hammingcount=0
        for i in range(len(leftlist)):
            
            if similaritycheck(selected,leftlist[i]) > 0:
                
                hammingcount+=1
        if hammingcount == len(leftlist): 
            if similaritycheck(selected,Complimentmap(selected)) > similaritycounter:
                if selected not in leftlist and Complimentmap(selected) not in leftlist:
                    leftlist.append(selected)
                    leftlist.append(Complimentmap(selected))
                    countleft+=1
    rightlist=[]
    countright=0
    while countright<motiveamounts:
        selected=generate_balanced_string(runthroughsizehalf,absdiff)
        hammingcount=0
        for i in range(len(rightlist)):
            if similaritycheck(selected,rightlist[i]) > 0:
                hammingcount+=1
        if hammingcount == len(rightlist): 
            if similaritycheck(selected,Complimentmap(selected)) > similaritycounter:
                if selected not in rightlist and Complimentmap(selected) not in rightlist and selected not in leftlist and Complimentmap(selected) not in leftlist:
                    rightlist.append(selected)
                    rightlist.append(Complimentmap(selected))
                    countright+=1

        finalset= [''.join(p) for p in itertools.product(leftlist, rightlist)]

        with open(f'./tests/testlibraries/lib_simple_{oligosize}bp_{amountofoligos}.csv', 'w') as file:
            for i in range(len(finalset)):
                file.write(f"{finalset[i]}\n")


