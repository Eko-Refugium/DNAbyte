import csv
import scipy.special
import math as m
import itertools as it
import random
import traceback
from operator import itemgetter
from pyxdameraulevenshtein import damerau_levenshtein_distance

from dnabyte.encode import Encode
from dnabyte.auxiliary import complementmap
from dnabyte.oligo import create_positional_libraries
from dnabyte.auxiliary import undoreedsolomon, undoltcode, check_first_n_entries_are_zero, binarize_and_combine_first_n_elements, remove_complements, find_closest_string
from dnabyte.oligo import create_positional_libraries

class DecodeLinearBinom(Encode):
    """
    This class provides binomial encoding for linear assembly based on the specified parameters.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def binary_list_to_decimal(self, binary_list):
        """
        Converts a list of binary digits to a decimal number.
        """
        # Convert the list of binary digits to a binary string
        binary_str = ''.join(str(digit) for digit in binary_list)
        # Convert the binary string to a decimal number
        decimal = int(binary_str, 2)
        return decimal
    

    def remove_complements(self,dna_list):
        """
        Removes complementary DNA sequences from a list, keeping only unique sequences.
        """
        unique_sequences = set()
        for seq in dna_list:
            comp_seq = complementmap(seq)
            if seq not in unique_sequences and comp_seq not in unique_sequences:
                unique_sequences.add(seq)
        
        return list(unique_sequences)

    def decode_binary_codewords(self, data, params):

        max_codeword_length = len(params.library.leftmotives) // 2 + len(params.library.rightmotives) // 2 - 1
        positionallibs = create_positional_libraries(params.library, len(params.library.leftmotives) + len(params.library.rightmotives) - 1)
        onepositionallib = positionallibs[0]
        oligo_length = len(onepositionallib[0])
        leftmotifofpos = []
        rightmotifofpos = []
        
        for j in range(len(positionallibs[0])):
            leftmotifofpos.append(positionallibs[0][j][:oligo_length//2])
            rightmotifofpos.append(positionallibs[0][j][oligo_length//2:])
        leftmotifofpos = remove_complements(leftmotifofpos)
        rightmotifofpos = remove_complements(rightmotifofpos)
        leftmotifofpos = list(sorted(set(leftmotifofpos)))
        rightmotifofpos = list(sorted(set(rightmotifofpos)))
        newpositionallibs =[''.join(combo) for combo in it.product(leftmotifofpos, rightmotifofpos)]

        listofzeros = [0] * len(newpositionallibs)
        for j in range(params.sigma_amount):
            listofzeros[j] = 1
        counterofwhereblock = 0
        
        while listofzeros != [0] * (len(listofzeros) - listofzeros.count(1)) + [1] * listofzeros.count(1):
            dnaswithoutcomplimetns = []
            for dnasinstep in range(len(listofzeros)):
                if listofzeros[dnasinstep] == 1:
                    dnaswithoutcomplimetns.append(newpositionallibs[dnasinstep])
            uniquemotivesleft = []
            uniquemotivesright = []
            for makinguniquecount in range(len(dnaswithoutcomplimetns)):
                uniquemotivesleft.append(dnaswithoutcomplimetns[makinguniquecount][:oligo_length//2])
                uniquemotivesright.append(dnaswithoutcomplimetns[makinguniquecount][oligo_length//2:])
            uniquemotivesleft = list(set(uniquemotivesleft))
            uniquemotivesright = list(set(uniquemotivesright))
            uniquecountleft = len(uniquemotivesleft)
            uniquecountright = len(uniquemotivesright)
            uniquecountoligos = uniquecountleft + uniquecountright
            counterofwhereblock = counterofwhereblock + 2 ** uniquecountoligos
            for chekforonetomove in range(len(listofzeros) - 1, -1, -1):
                if listofzeros[chekforonetomove] == 0 and listofzeros[chekforonetomove - 1] == 1:
                    l = chekforonetomove -1
                    listofzeros[l] = 0
                    listofzeros[l + 1] = 1
                    ones_to_move = listofzeros[l+2:]
                    listofzeros[l + 2:] = [0] * len(ones_to_move)
                    listofzeros[l + 2:] = [1] * ones_to_move.count(1) + [0] * ones_to_move.count(0)
                    break
        
        if params.inner_error_correction == 'ltcode':
            message_length = max_codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1 - params.index_carry_length - params.ltcode_header
        else:
            message_length = max_codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1

        if params.outer_error_correction == 'reedsolomon':
            positionallibs = create_positional_libraries(params.library, len(params.library.leftmotives) + len(params.library.rightmotives) - 1)
            onepositionallib = positionallibs[0]
            fixingeverything = []
            howmanybitsin = m.floor(m.log2(counterofwhereblock))
            keepgoing = True
            while keepgoing:
                if howmanybitsin % 8 == 0:
                    keepgoing = False
                else:
                    howmanybitsin -= 1
            while message_length * howmanybitsin//8 > 255:
                message_length -= 1

            tempdata = []
            for i, codewords in enumerate(data):
                codewords = [x for x in codewords if len(bin(x)[2:].zfill(howmanybitsin)) == howmanybitsin]
                tempdata.append(codewords)
            tempdata = [x for x in tempdata if len(x) > 0]
            data = tempdata

            params.reed_solo_percentage = params.reed_solo_percentage
            reedsoloinfo = m.floor(message_length * params.reed_solo_percentage)
            reedsolocorr = message_length - reedsoloinfo
            if params.inner_error_correction == 'ltcode':
                undoneltcode, checkforvalidlt = undoltcode(data,params.index_carry_length, params.ltcode_header, howmanybitsin)
            else:
                undoneltcode = data
            
            decodeddata, checkforvalidreed = undoreedsolomon(undoneltcode, reedsolocorr, howmanybitsin)
            for codewords in decodeddata:
                
                if check_first_n_entries_are_zero(codewords, params.codeword_maxlength_positions):
                    codewords=codewords[params.codeword_maxlength_positions:]
                else:
                    howmanydoihavetopop = binarize_and_combine_first_n_elements(codewords,params.codeword_maxlength_positions, howmanybitsin)
                    codewords = codewords[: -howmanydoihavetopop]
                    codewords = codewords[params.codeword_maxlength_positions:]
                
                if check_first_n_entries_are_zero(codewords, 1):
                    codewords = codewords[1:]
                    binaries = [format(number, '0{}b'.format(howmanybitsin)) for number in codewords]
                else:
                    howmanydoihavetopop = binarize_and_combine_first_n_elements(codewords,1, howmanybitsin)
                    binaries = [format(number, '0{}b'.format(howmanybitsin)) for number in codewords]
                    binaries[-1] = binaries[-1][(howmanybitsin-howmanydoihavetopop):]
                    binaries = binaries[1:]
                
                fixingeverything.append(binaries)
            
            flattened_list = [item for sublist in fixingeverything for item in sublist]
            combined_string = ''.join(flattened_list)
            if params.inner_error_correction == 'ltcode':
                if checkforvalidlt and checkforvalidreed:
                    checkall = True
                else:
                    checkall = False
            else:
                if checkforvalidreed:
                    checkall = True
                else:   
                    checkall = False
            return combined_string, checkall
        else:
            positionallibs = create_positional_libraries(params.library, len(params.library.leftmotives)+len(params.library.rightmotives) - 1)
            onepositionallib = positionallibs[0]
            howmanybitsin = m.floor(m.log2(counterofwhereblock))
            fixingeverything = []
            if params.inner_error_correction == 'ltcode':
                undoneltcode, checkforvalidlt = undoltcode(data, params.index_carry_length, params.ltcode_header, howmanybitsin)
            else:
                undoneltcode = data
            for codewords in undoneltcode:
                
                if check_first_n_entries_are_zero(codewords, params.codeword_maxlength_positions):
                    codewords = codewords[params.codeword_maxlength_positions:]
                else:
                    howmanydoihavetopop = binarize_and_combine_first_n_elements(codewords, params.codeword_maxlength_positions, howmanybitsin)
                    codewords = codewords[:-howmanydoihavetopop]
                    codewords = codewords[params.codeword_maxlength_positions:]
                if check_first_n_entries_are_zero(codewords, 1):
                    codewords = codewords[1:]
                    binaries = [format(number, '0{}b'.format(howmanybitsin)) for number in codewords]
                else:
                    howmanydoihavetopop=binarize_and_combine_first_n_elements(codewords,1, howmanybitsin)
                    binaries = [format(number, '0{}b'.format(howmanybitsin)) for number in codewords]
                    binaries[-1] = binaries[-1][(howmanybitsin - howmanydoihavetopop):]
                    binaries = binaries[1:]
                
                fixingeverything.append(binaries)
            
            flattened_list = [item for sublist in fixingeverything for item in sublist]
            combined_string = ''.join(flattened_list)
            if params.inner_error_correction == 'ltcode':
                if checkforvalidlt:
                    checkall = True
                else:
                    checkall = False
            else:
                checkall = True
            return combined_string, checkall


    def decode_linear_binom(self, data):
        """
        Decodes the provided data using binomial encoding.
        """
        try:
            DNAs = self.params.library.library
            DNAs.sort()
            oligo_length = len(DNAs[0])

            leftmotives = self.params.library.rightmotives
            rightmotives = self.params.library.leftmotives
            
            leftmotif = list(sorted(set(leftmotives)))
            rightmotif = list(sorted(set(rightmotives)))
            
            positionallibs = create_positional_libraries(self.params.library,  len(leftmotif) + len(rightmotif) - 1)
            for i in range(len(leftmotives)):
                leftmotives.append(complementmap(leftmotives[i]))
            for i in range(len(rightmotives)):
                rightmotives.append(complementmap(rightmotives[i]))
            
            binariesinchuncks = []
            for codeword in range(len(data.data)):
                binary = []
                counterofwhereblock = 0
                for i in range(len(data.data[codeword])):
                    leftmotifofpos = []
                    rightmotifofpos = []
                    allmotivesofsublist = []
                    comparemotives = []
                    sorted(positionallibs[i + self.params.dna_barcode_length])
                    for j in range(len(positionallibs[i + self.params.dna_barcode_length])):
                        leftmotifofpos.append(positionallibs[i + self.params.dna_barcode_length][j][:oligo_length // 2])
                        rightmotifofpos.append(positionallibs[i + self.params.dna_barcode_length][j][oligo_length // 2:])

                    sorted(leftmotifofpos)
                    sorted(rightmotifofpos)
                        
                    leftmotifofpos = self.remove_complements(leftmotifofpos)
                    rightmotifofpos = self.remove_complements(rightmotifofpos)

                    leftmotifofpos = list(sorted(set(leftmotifofpos)))
                    rightmotifofpos = list(sorted(set(rightmotifofpos)))

                    allmotivesofsublist = list(set(allmotivesofsublist))

                    for k in range(len(leftmotifofpos)):
                        allmotivesofsublist.append(leftmotifofpos[k])

                    for k in range(len(rightmotifofpos)):
                        allmotivesofsublist.append(rightmotifofpos[k])

                    newpositionallibs =[''.join(combo) for combo in it.product(leftmotifofpos, rightmotifofpos)]
                                    
                    uniquemotivesleft = []
                    uniquemotivesright = []
                    uniquemotivesalphabetleft = []
                    uniquemotivesalphabetright = []
                    alluniquemotives = []
                    alluniquemotivesalphabet = []

                    for makinguniquecount in range(len(data.data[codeword][i])):
                        uniquemotivesleft.append(data.data[codeword][i][makinguniquecount][:oligo_length // 2])
                        uniquemotivesalphabetleft.append(data.data[codeword][i][makinguniquecount][:oligo_length // 2])
                        uniquemotivesright.append(data.data[codeword][i][makinguniquecount][oligo_length // 2:])
                        uniquemotivesalphabetright.append(data.data[codeword][i][makinguniquecount][oligo_length // 2:])

                    for k in range(len(leftmotifofpos)):
                        if leftmotifofpos[k] in uniquemotivesleft or complementmap(leftmotifofpos[k]) in uniquemotivesleft:
                            comparemotives.append(leftmotifofpos[k])

                    for k in range(len(rightmotifofpos)):
                        if rightmotifofpos[k] in uniquemotivesright or complementmap(rightmotifofpos[k]) in uniquemotivesright:
                            comparemotives.append(rightmotifofpos[k])

                    for k in range(len(uniquemotivesalphabetleft)):
                        if uniquemotivesalphabetleft[k] > complementmap(uniquemotivesalphabetleft[k]):
                            uniquemotivesalphabetleft[k] = complementmap(uniquemotivesalphabetleft[k])

                    for k in range(len(uniquemotivesalphabetright)):
                        if uniquemotivesalphabetright[k] > complementmap(uniquemotivesalphabetright[k]):
                            uniquemotivesalphabetright[k] = complementmap(uniquemotivesalphabetright[k])
                    
                    uniquemotivesleft = list(set(uniquemotivesleft))
                    uniquemotivesright = list(set(uniquemotivesright))
                    uniquemotivesalphabetleft = list(set(uniquemotivesalphabetleft))
                    uniquemotivesalphabetright = list(set(uniquemotivesalphabetright))

                    complementuniquemotivesleft = []
                    complementuniquemotivesright = []

                    for k in range(len(uniquemotivesleft)):
                        if complementmap(uniquemotivesleft[k]) not in complementuniquemotivesleft and uniquemotivesleft[k] not in complementuniquemotivesleft:
                            complementuniquemotivesleft.append(uniquemotivesleft[k])

                    for k in range(len(uniquemotivesright)):
                        if complementmap(uniquemotivesright[k]) not in complementuniquemotivesright and uniquemotivesright[k] not in complementuniquemotivesright:
                            complementuniquemotivesright.append(uniquemotivesright[k])

                    uniquecountleft = len(complementuniquemotivesleft)
                    uniquecountright = len(complementuniquemotivesright)

                    for k in range(len(complementuniquemotivesleft)):
                        alluniquemotives.append(complementuniquemotivesleft[k])

                    for k in range(len(complementuniquemotivesright)):
                        alluniquemotives.append(complementuniquemotivesright[k])

                    for k in range(len(uniquemotivesalphabetleft)):
                        alluniquemotivesalphabet.append(uniquemotivesalphabetleft[k])

                    for k in range(len(uniquemotivesalphabetright)):
                        alluniquemotivesalphabet.append(uniquemotivesalphabetright[k])
                    
                    uniquecountoligos = uniquecountleft+uniquecountright
                    sortedlsit = []

                    for k in range(len(comparemotives)):
                        if comparemotives[k] in alluniquemotives:
                            sortedlsit.append(comparemotives[k])
                        elif complementmap(comparemotives[k]) in alluniquemotives:
                            sortedlsit.append(complementmap(comparemotives[k]))
                        
                    thecompliments = [0]*uniquecountoligos

                    for uniquecount in range(uniquecountoligos):

                        if sortedlsit[uniquecount] not in allmotivesofsublist:
                            thecompliments[uniquecount] = 1
                
                    extradecimal = self.binary_list_to_decimal(thecompliments)
                    
                    for turningcompliments in range(len(sortedlsit)):
                        if turningcompliments < uniquecountleft and thecompliments[turningcompliments] == 1:
                            for checkingcompliments in range(len(data.data[codeword][i])):
                                if data.data[codeword][i][checkingcompliments][:oligo_length // 2] == sortedlsit[turningcompliments]:
                                    data.data[codeword][i][checkingcompliments] = complementmap(data.data[codeword][i][checkingcompliments][:oligo_length//2])+data.data[codeword][i][checkingcompliments][oligo_length//2:]
                        if turningcompliments >= uniquecountleft and thecompliments[turningcompliments] == 1:
                            for checkingcompliments in range(len(data.data[codeword][i])):
                                if data.data[codeword][i][checkingcompliments][oligo_length // 2:] == sortedlsit[turningcompliments]:
                                    data.data[codeword][i][checkingcompliments] = data.data[codeword][i][checkingcompliments][:oligo_length//2]+complementmap(data.data[codeword][i][checkingcompliments][oligo_length//2:]) 
                    makingtheoneslisz = [0] * len(newpositionallibs)

                    for l in range(len(data.data[codeword][i])):

                        if data.data[codeword][i][l] not in newpositionallibs:
                            runningloop = True
                            tempcodewordpsrt = self.find_closest_string(data.data[codeword][i][l], newpositionallibs)
                            if tempcodewordpsrt in newpositionallibs and tempcodewordpsrt not in data.data[codeword][i]:
                                data.data[codeword][i][l]=tempcodewordpsrt
                            else:
                                while runningloop:
                                    randomoligo = random.randint(0, len(newpositionallibs) - 1)
                                    tempcodewordpsrt = newpositionallibs[randomoligo]
                                    if tempcodewordpsrt in newpositionallibs and tempcodewordpsrt not in data.data[codeword][i]:
                                        data.data[codeword][i][l] = tempcodewordpsrt
                                        runningloop = False
                                    
                        makingtheoneslisz[newpositionallibs.index(data.data[codeword][i][l])] = 1
                    if sum(makingtheoneslisz) != self.params.sigma_amount:
                        while sum(makingtheoneslisz) != self.params.sigma_amount:
                            randomzeroposition = random.randint(0,len(makingtheoneslisz)-1)
                            if makingtheoneslisz[randomzeroposition] == 0:
                                makingtheoneslisz[randomzeroposition] = 1
                                        
                    listofzeros = [0] * len(newpositionallibs)
                    
                    for j in range(self.params.sigma_amount):
                        listofzeros[j] = 1
                    looper = True
                    
                    counterofwhereblock = 0
                    while looper:
                        dnaswithoutcomplimetns = []
                        for dnasinstep in range(len(listofzeros)):
                            if listofzeros[dnasinstep] == 1:
                                dnaswithoutcomplimetns.append(newpositionallibs[dnasinstep])
                        uniquemotivesleft = []
                        uniquemotivesright = []
                        for makinguniquecount in range(len(dnaswithoutcomplimetns)):
                            uniquemotivesleft.append(dnaswithoutcomplimetns[makinguniquecount][:oligo_length // 2])
                            uniquemotivesright.append(dnaswithoutcomplimetns[makinguniquecount][oligo_length // 2:])
                        
                        uniquemotivesleft = list(sorted(set(uniquemotivesleft)))
                        uniquemotivesright = list(sorted(set(uniquemotivesright)))

                        uniquecountleft = len(uniquemotivesleft)
                        uniquecountright = len(uniquemotivesright)

                        uniquecountoligos = uniquecountleft + uniquecountright
                        lastblockborder = counterofwhereblock
                        counterofwhereblock = counterofwhereblock + 2 ** uniquecountoligos
                        if listofzeros == makingtheoneslisz:
                            looper = False
                        else: 
                            for chekforonetomove in range(len(listofzeros)-1, -1, -1):
                                if listofzeros[chekforonetomove] == 0 and listofzeros[chekforonetomove - 1] == 1:
                                    l = chekforonetomove -1
                                    listofzeros[l] = 0
                                    listofzeros[l + 1] = 1
                                    ones_to_move = listofzeros[l+2:]
                                    listofzeros[l + 2:] = [0] * len(ones_to_move)
                                    listofzeros[l + 2:] = [1] * ones_to_move.count(1) + [0] * ones_to_move.count(0)
                                    break
                    
                    sumer = lastblockborder + extradecimal 
                    binary.append(int(sumer))
                binariesinchuncks.append(binary)
            
            decoded_binary, valid = self.decode_binary_codewords(binariesinchuncks, self.params)
    
            info = {
                "number_of_codewords": len(decoded_binary),
                "data_length": len(data.data)
            }
        
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error during decoding: {e}")
                self.logger.error(traceback.format_exc())

            decoded_binary = None
            valid = False
            info = None

        return decoded_binary, valid, info

        