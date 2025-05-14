import csv
import scipy.special
import math as m
import itertools as it
import random
from dnabyte.auxiliary import complementmap


def decimal_to_binary(n): 
    return bin(n).replace("0b", "") 

def transpose_lists(list_of_lists):
    # Use the zip function with argument unpacking to transpose the list of lists
    transposed = [list(group) for group in zip(*list_of_lists)]
    return transposed

inputbinary = "1111111111111111101101010010101001010010100101010010100101010010101010010111110101010001010101111111101111110010101001010100101010010101001010110100101010101010101011101010010100010111011111010101011111111101111110111101001010100101010101001010100000101010111111111110101010110000001010101111111111"

#importing library

###rongwells=[]
DNAs=[]
with open('20bp_Lib.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        
        if line_count == 0:
            line_count+=1
        else:
           # wrongwells.append(row[0])
            DNAs.append(row[2])
            line_count+=1


#finding all left and right motives of the dnas

leftmotives = []
rightmotives = []
for i in range(len(DNAs)):
    lengthofDNA = len(DNAs[i])
    leftpart = DNAs[i][:lengthofDNA//2]
    rightpart = DNAs[i][lengthofDNA//2:]
    if leftpart not in leftmotives and complementmap(leftpart) not in leftmotives:
        leftmotives.append(leftpart)
    if rightpart not in rightmotives and complementmap(rightpart) not in rightmotives:
        rightmotives.append(rightpart)
leftmotives.sort()
rightmotives.sort()

# parameters
lengthofcodeword = 7
#Here i need to fillter the dnas
correctDNAs = DNAs
correctDNAs.sort()
sigmaamount = 5
coppynumbers = 100

positionallibrary = []
for remover in range(len(correctDNAs)):
    #check for rightmotives?
    if leftmotives[0] not in correctDNAs[remover] and complementmap(leftmotives[0]) not in correctDNAs[remover] and rightmotives[0] not in correctDNAs[remover] and complementmap(rightmotives[0]) not in correctDNAs[remover]:
        positionallibrary.append(correctDNAs[remover])
positionallibrary.sort() 

binomcoeff = m.floor(m.log2(scipy.special.binom(len(positionallibrary), sigmaamount)))
if binomcoeff%2 != 0:
    binomcoeff = binomcoeff-1
    
howmanybitsin = binomcoeff

# making binary blocks
binaryblocks = []
for i in range(0, len(inputbinary), howmanybitsin):
    binaryblocks.append(inputbinary[i:i+howmanybitsin])
safelengthofthelastblock = len(binaryblocks[-1])
binaryblocks[-1]=binaryblocks[-1].zfill(howmanybitsin) #filling the last block with zeros

# making the binary codewords
binarycodewords = []
for i in range(0, len(binaryblocks), lengthofcodeword-1):
    counter = 0
    onecodeowrd = []
    while counter < lengthofcodeword:
        if i+counter >= len(binaryblocks):
            break
        onecodeowrd.append(binaryblocks[i+counter])
        counter+=1
    binarycodewords.append(onecodeowrd)
for i in range(len(binarycodewords)):
    binarycodewords[i].insert(0, str(decimal_to_binary(i)).zfill(howmanybitsin))

#making the using library
DNAcombinatorialcode = []
decimalchecker = []
for overall in range(len(binarycodewords)):
    tempitemp = []
    for i in range(len(binarycodewords[overall])):
        positionallibrary = []
        for remover in range(len(correctDNAs)):
            #check for rightmotives?
            if leftmotives[i] not in correctDNAs[remover] and complementmap(leftmotives[i]) not in correctDNAs[remover] and rightmotives[i] not in correctDNAs[remover] and complementmap(rightmotives[i]) not in correctDNAs[remover]:
                positionallibrary.append(correctDNAs[remover])
            #Check if sticky ends are correct
        positionallibrary.sort() 


        #making the translation table
        # this takes ages
        # lookuptable= list(it.combinations(positionallibrary, sigmaamount))
            
        combinatorialcodeword = []
        combinatorialcodewordoneolig = []
        listofzeros = [0] * len(positionallibrary)
        decimal= int(binarycodewords[overall][i],2)
        decimalchecker.append(decimal)
        comparevalue=0
        checker=0
        for k in range(sigmaamount):
            firstnumberchecker = False
            
            while not firstnumberchecker:
                comparevalue = comparevalue + scipy.special.binom(len(positionallibrary)-checker-1, sigmaamount-k-1)
                if decimal<comparevalue:
                    firstnumberchecker = True
                    listofzeros[checker]=1
                    comparevalue = comparevalue - scipy.special.binom(len(positionallibrary)-checker-1, sigmaamount-k-1)
                checker+=1
        combinatorialcodewordoneolig.append(listofzeros)
        DNAcombinatorialcodetemp = []
        for j in range(len(combinatorialcodewordoneolig)):
            DNAcombinatorialcodetemptemp = []
            for k in range(len(combinatorialcodewordoneolig[j])):
                if combinatorialcodewordoneolig[j][k] == 1:
                    DNAcombinatorialcodetemptemp.append(positionallibrary[k])
            #DNAcombinatorialcodetemp.append(DNAcombinatorialcodetemptemp)
        tempitemp.append(DNAcombinatorialcodetemptemp)
    DNAcombinatorialcode.append(tempitemp)
    
oligotrios = []
for i in range(len(DNAcombinatorialcode)):
    oligotriostemp0 = []
    for j in range(len(DNAcombinatorialcode[i])):
        oligotriostemp1 = []
        for k in range(len(DNAcombinatorialcode[i][j])):
            oligotriostemp2 = []
        
            if j%2==0:
                oligotriostemp2.append(leftmotives[j//2]+complementmap(DNAcombinatorialcode[i][j][k][len(correctDNAs[0])//2:]))
                oligotriostemp2.append(DNAcombinatorialcode[i][j][k])
                oligotriostemp2.append(complementmap(DNAcombinatorialcode[i][j][k][:len(correctDNAs[0])//2])+rightmotives[j//2])
            else:
                oligotriostemp2.append(complementmap(DNAcombinatorialcode[i][j][k][:len(correctDNAs[0])//2])+complementmap(rightmotives[j//2]))
                oligotriostemp2.append((DNAcombinatorialcode[i][j][k]))
                oligotriostemp2.append(complementmap(leftmotives[j//2])+complementmap(DNAcombinatorialcode[i][j][k][len(correctDNAs[0])//2:]))
            
            oligotriostemp1.append(oligotriostemp2)
        oligotriostemp0.append(oligotriostemp1)
    oligotrios.append(oligotriostemp0)


#making the mixtures
mixtures = []
for j in range(len(oligotrios)):
    tempsplit = []
    for i in range(coppynumbers):
        tempsmixutres = []
        for k in range(len(oligotrios[j])):
            randomcoise = random.choice(oligotrios[j][k])
            for l in range(len(randomcoise)):
                tempsmixutres.append(randomcoise[l])
        tempsplit.append(tempsmixutres)
    mixtures.append(tempsplit)

fiveto3strings = []
for i in range(len(mixtures)):
    tempsplitfull = []
    for k in range(len(mixtures[i])):
        
        fiveto3string = ''
        for j in range(len(mixtures[i][k])):
            if j%2==0:
                fiveto3string+=mixtures[i][k][j]
        tempsplitfull.append(fiveto3string)
    fiveto3strings.append(tempsplitfull)

#Decoding

#extracting data

oligolength = len(correctDNAs[0])
informationoligos = []
for i in range(len(fiveto3strings)):
    tempinformationoligos = []
    for j in range(len(fiveto3strings[i])):
        infomationishere=[]
        for k in range(oligolength//2,len(fiveto3strings[i][j]),oligolength//2):
            if k%3!=0:
                ifostringbuilder = ''
                for l in range(oligolength//2):
                    ifostringbuilder+=fiveto3strings[i][j][k+l]
                infomationishere.append(ifostringbuilder)
        
        combined_pairs = []
        for fixer in range(0, len(infomationishere), 2):
            if fixer+1 < len(infomationishere):
                combined = infomationishere[fixer] + infomationishere[fixer+1]
                combined_pairs.append(combined)
            else:
                combined_pairs.append(infomationishere[fixer])
        for flipper in range(len(combined_pairs)):
            if flipper%2==0:
                combined_pairs[flipper] = complementmap(combined_pairs[flipper])
        tempinformationoligos.append(combined_pairs)
    informationoligos.append(tempinformationoligos)

#finding the 5 oligos

posstionalshift = []
for i in range(len(informationoligos)):
    posstionalshift.append(transpose_lists(informationoligos[i]))
for i in range(len(posstionalshift)):
    for j in range(len(posstionalshift[i])):
        posstionalshift[i][j] = list(set(posstionalshift[i][j]))

# for i in range(len(binarycodewords)):
#     positionallibrary = []
#     for remover in range(len(correctDNAs)):
#         #check for rightmotives?
#         if leftmotives[i] not in correctDNAs[remover] and complementmap(leftmotives[i]) not in correctDNAs[remover] and rightmotives[i] not in correctDNAs[remover] and complementmap(rightmotives[i]) not in correctDNAs[remover]:
#             positionallibrary.append(correctDNAs[remover])
#     positionallibrary.sort() 

#     listofrecoveries = []
    
#     for j in range(len(posstionalshift)):
#         temprecovery = []
#         makingtheoneslisz= [0] * len(positionallibrary)
#         for l in range(len(posstionalshift[j][i])):
#             for n in range(len(positionallibrary)):
#                 if posstionalshift[i][j][l] == positionallibrary[n]:
#                     makingtheoneslisz[n] = 1
#         temprecovery.append(makingtheoneslisz)
#         listofrecoveries.append(temprecovery)


listofrecoveries = []
for i in range(len(posstionalshift)):
    tamprecovery = []
    
    for j in range(len(posstionalshift[i])):

        makingtheoneslisz= [0] * len(positionallibrary)
        positionallibrary = []
        for remover in range(len(correctDNAs)):
            #check for rightmotives?
            if leftmotives[j] not in correctDNAs[remover] and complementmap(leftmotives[j]) not in correctDNAs[remover] and rightmotives[j] not in correctDNAs[remover] and complementmap(rightmotives[j]) not in correctDNAs[remover]:
                positionallibrary.append(correctDNAs[remover])
        positionallibrary.sort() 

        for l in range(len(posstionalshift[i][j])):
            for n in range(len(positionallibrary)):
                if posstionalshift[i][j][l] == positionallibrary[n]:
                    makingtheoneslisz[n] = 1
        counterforsum =0
        sumer = 0
        onescounter=0
        for camlulator in range(len(makingtheoneslisz)):
            if onescounter != sigmaamount:
                
                if makingtheoneslisz[camlulator] == 0:
                    sumer+=scipy.special.binom(len(positionallibrary)-camlulator-1, sigmaamount-onescounter-1)
                    
                else:
                    onescounter+=1
        counterforsum+=1     
            
        tamprecovery.append(sumer)
    listofrecoveries.append(tamprecovery)


sorted_list = sorted(listofrecoveries, key=lambda x: x[0])
modified_listofrecoveries = [sublist[1:] for sublist in listofrecoveries]
flattened_list = [item for sublist in modified_listofrecoveries for item in sublist]
for i in range(len(flattened_list)):
    flattened_list[i] = int(flattened_list[i])
binaries = [format(number, '0{}b'.format(howmanybitsin)) for number in flattened_list]
combined_string = ''.join(binaries)

