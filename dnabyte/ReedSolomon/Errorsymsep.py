import random
import math as m
import numpy as np

def GCerr(ErrorFuncIns, ErrorFuncDel, ErrorFuncSub, InputDNA): #GCerror calcualtion function
    DNABases=np.array(4)
    DNABases=['A','G','C','T']
    DNABasesn=np.array(4)
    DNABasesn=range(255)

    unique, counts = np.unique(InputDNA, return_counts=True)
    #Gamount=InputDNA.count(1)
    #Camount=InputDNA.count(2)
    gind=np.where(unique==1)[0]
    cind=np.where(unique==2)[0]
    GCcont=(counts[gind]+counts[cind])/len(InputDNA)
    ProbIns=ErrorFuncIns(GCcont) #Very simplified rn one would need ro do it iteratively to find the parts of the DNA that actually contains a extrem GC count
    ProbDel=ErrorFuncDel(GCcont)
    ProbSub=ErrorFuncSub(GCcont)   
    OutputDNA=InputDNA
    endofloop=len(InputDNA)
    
    for i in range(len(InputDNA)): #exicuting the acctual errors
        nIn = random.random()
        nDel = random.random()
        nSub = random.random()

        if ProbIns>nIn:
            random_index = random.randint(0, len(DNABases) - 1)
            random_element = list(DNABases[random_index])
            OutputDNA=np.insert(OutputDNA,i,random_element)
            
        
        if ProbDel>nDel:
            OutputDNA=np.delete(OutputDNA,i)
            endofloop=endofloop-1

        if (1-(1-ProbSub)**4)>nSub:
            what=InputDNA[i]
            random_index = random.randint(0, len(DNABasesn)-1)
            random_element = DNABasesn[random_index]
            OutputDNA[i]=random_element
            DNABasesn=range(255)

    return OutputDNA