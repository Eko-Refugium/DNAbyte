import numpy as np

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

def fountaindecode(data,listlen,sizeofmesseage):
    sorted_list = sorted(data, key=lambda x: x.degree)
    finallist,restlist=initilation(sorted_list)
    finfinlisttemp=reprodice(finallist,restlist,listlen)
    finfinlist=[]
    for i in finfinlisttemp:
        dropbin=np.base_repr(int(i.data))
        dropbin= '0' * (sizeofmesseage - len(dropbin)) + dropbin
        finfinlist.append(translate_binary(dropbin))
    finfinlist=sorted(list(set(finfinlist)))
    return finfinlist