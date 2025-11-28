import galois
import random
import rsKaya
import base_conversion_kit as bc
import Errorsymsep as err
import numpy as np

def base_convert(i, b):
    result = []
    while i > 0:
            result.insert(0, i % b)
            i = i // b
    
    while len(result)!=4:
        result.insert(0,0)
    
    return result

def solve( A,  B, stringlen):
  
    count = 0
 
    # since, the numbers are less than 2^31
    # run the loop from '0' to '31' only
    for i in range(0,stringlen):
 
        # right shift both the numbers by 'i' and
        # check if the bit at the 0th position is different
        if ((( A >>  i) & 1) != (( B >>  i) & 1)): 
             count=count+1
    return count

dir(galois)
k=223
octallen=k*4
octnum=np.zeros(int(octallen))
for i in range(octallen):
    octnum[i]=int(random.randint(0, 3))


coder=rsKaya.RScoderKaya(255,223)
temphold=[]
for i in range(int(len(octnum)/4)):
    temstringgrp=''
    for eight in range(4):
        temstringgrp=temstringgrp+str(int(octnum[i*4+eight]))
    temphold.append(int(temstringgrp))  
encoded_b10=[]
for i in range(len(temphold)):
    encoded_b10.append(int(bc.convert_base(int(temphold[i]),4,10)))
c=coder.generate(encoded_b10)


def fInGC(x):
    return 0

def fDelGC(x):
    return 0

def fSubGC(x):
    return 0.1


c[5]=0
c[7]=2
c[36]=0
c[22]=0


countererr=0
posoferrre=[]
for posi in range(len(encoded_b10)):
     if c[posi]!=encoded_b10[posi]:
        countererr=countererr +1
        posoferrre.append(posi)


f=coder.correct(c)

counter=0
posoferr=[]
for posi in range(len(f)):
     if f[posi]!=octnum[posi]:
        counter=counter +1
        posoferr.append(posi)
