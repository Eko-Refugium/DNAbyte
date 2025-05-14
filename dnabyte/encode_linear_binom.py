import scipy.special
import math as m
import itertools as it
import random
import traceback

from dnabyte.oligo import create_positional_libraries
from dnabyte.auxiliary import create_counter_list, complementmap
from dnabyte.encode import Encode
from dnabyte.auxiliary import decimal_to_binary, MakeReedSolomonCode, makeltcode, remove_complements
from dnabyte.oligo import create_positional_libraries

class EncodeLinearBinom(Encode):
    """
    This class provides binomial encoding for linear assembly based on the specified parameters.
    """
    def __init__(self, params, logger=None):
        self.params = params
        self.logger = logger

    def combine_elements_randomized(self, list_of_lists):
        """
        Takes a list of lists and generates a list of lists where each new list contains one element from each sublist,
        ensuring that each element from the original sublists is used exactly once, in a randomized order.

        :param list_of_lists: The list of lists to combine.
        :return: A list of lists with combined elements in a randomized order.
        """
        # Ensure all sublists have the same length
        if not all(len(sublist) == len(list_of_lists[0]) for sublist in list_of_lists):
            raise ValueError("All sublists must have the same length")

        # Number of elements in each sublist
        n = len(list_of_lists[0])

        # Generate a list of indices and shuffle it
        indices = list(range(n))
        random.shuffle(indices)

        # Generate the combined lists using the shuffled indices
        combined_lists = []
        for i in indices:
            combined_list = [sublist[i] for sublist in list_of_lists]
            combined_lists.append(combined_list)

        return combined_lists
    

    def remove_complements(self, dna_list):
        unique_sequences = set()
        for seq in dna_list:
            comp_seq = complementmap(seq)
            if seq not in unique_sequences and comp_seq not in unique_sequences:
                unique_sequences.add(seq)
        
        return list(unique_sequences)

    def create_binary_codewords(self, data, params):

        oligo_length = len(params.library.library[0])

        # extract the left and right motifs from the positional library
        positionallibs = create_positional_libraries(params.library, len(params.library.leftmotives) + len(params.library.rightmotives) - 1)
        leftmotifofpos = []
        rightmotifofpos = []

        for j in range(len(positionallibs[0])):
            leftmotifofpos.append(positionallibs[0][j][:oligo_length // 2])
            rightmotifofpos.append(positionallibs[0][j][oligo_length // 2:])

        leftmotifofpos = remove_complements(leftmotifofpos)
        rightmotifofpos = remove_complements(rightmotifofpos)
        leftmotifofpos = list(set(leftmotifofpos))
        rightmotifofpos = list(set(rightmotifofpos))
        leftmotifofpos.sort()
        rightmotifofpos.sort()

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
            for chekforonetomove in range(len(listofzeros)-1, -1, -1):
                if listofzeros[chekforonetomove] == 0 and listofzeros[chekforonetomove - 1] == 1:
                    l = chekforonetomove - 1
                    listofzeros[l] = 0
                    listofzeros[l + 1] = 1
                    ones_to_move = listofzeros[l+2:]
                    listofzeros[l+2:] = [0] * len(ones_to_move)
                    listofzeros[l+2:] = [1] * ones_to_move.count(1) + [0] * ones_to_move.count(0)
                    break
        

        if params.outer_error_correction == 'reedsolomon':
            howmanybitsin = m.floor(m.log2(counterofwhereblock))
            keepgoing = True
            while keepgoing:
                if howmanybitsin % 8 == 0:
                    keepgoing = False
                else:
                    howmanybitsin-=1
        else:
            howmanybitsin = m.floor(m.log2(counterofwhereblock))
        

        # create the binary blocks
        binary_blocks = []
        for i in range(0, len(data), howmanybitsin):
            binary_blocks.append(data[i:i + howmanybitsin])
        number_of_zero_padding = len(binary_blocks[-1])
        binary_blocks[-1] = binary_blocks[-1].zfill(howmanybitsin)

        # TODO: rename the variable
        if number_of_zero_padding == howmanybitsin:
            number_of_zero_padding = 0


        max_codeword_length = len(params.library.leftmotives) + len(params.library.rightmotives) - 1

        if params.inner_error_correction == 'ltcode':
            message_length = max_codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1 - params.index_carry_length - params.ltcode_header
        else:
            message_length = max_codeword_length - params.dna_barcode_length - params.codeword_maxlength_positions - 1

        if params.outer_error_correction == 'reedsolomon':
            
            while message_length * howmanybitsin//8 > 255:
                message_length -= 1
            reedsoloinfo = m.floor(message_length * params.reed_solo_percentage)
            reedsolocorr = message_length - reedsoloinfo
            binarycodewords = []
            for i in range(0, len(binary_blocks), reedsoloinfo):
                counter = 0
                one_codeword = []
                while counter < reedsoloinfo:
                    if i + counter >= len(binary_blocks):
                        break
                    one_codeword.append(binary_blocks[i+counter])
                    counter += 1
                
            
                one_codeword.insert(0, str(decimal_to_binary(0)).zfill(howmanybitsin))
                counterhowmanybytesadded = 0
                while len(one_codeword) <= reedsoloinfo:
                    one_codeword.append('0' * howmanybitsin)
                    counterhowmanybytesadded += 1

                s = str(decimal_to_binary(counterhowmanybytesadded)).zfill(howmanybitsin * params.codeword_maxlength_positions)
                bitsofaddedcounter=[s[i:i + howmanybitsin] for i in range(0, len(s), howmanybitsin)]
                for i in reversed(range(len(bitsofaddedcounter))):
                    one_codeword.insert(0, bitsofaddedcounter[i].zfill(howmanybitsin))

                binarycodewords.append(one_codeword)
                if len(binary_blocks)//(reedsoloinfo)!= 0:
                    ValueError('The last codeword is not the right length')



            binarycodewords[-1][len(bitsofaddedcounter)] = str(decimal_to_binary(number_of_zero_padding)).zfill(howmanybitsin)           
            codewordatgoodlength = MakeReedSolomonCode(binarycodewords, reedsolocorr, howmanybitsin)  
            uniformnameforoutput = codewordatgoodlength
        
            

        else:
            codewordlength = message_length
            binarycodewords = []
            for i in range(0, len(binary_blocks), codewordlength):
                counter = 0
                one_codeword = []
                while counter < codewordlength:
                    if i+counter >= len(binary_blocks):
                        break
                    one_codeword.append(binary_blocks[i+counter])
                    counter += 1
                
                one_codeword.insert(0, str(decimal_to_binary(0)).zfill(howmanybitsin))
                counterhowmanybytesadded = 0
                while len(one_codeword) <= codewordlength:
                    one_codeword.append('0' * howmanybitsin)
                    counterhowmanybytesadded += 1
                s=str(decimal_to_binary(counterhowmanybytesadded)).zfill(howmanybitsin*params.codeword_maxlength_positions)
                bitsofaddedcounter=[s[i:i+howmanybitsin] for i in range(0, len(s), howmanybitsin)]
                for i in reversed(range(len(bitsofaddedcounter))):
                    one_codeword.insert(0, bitsofaddedcounter[i].zfill(howmanybitsin))

                binarycodewords.append(one_codeword)
                if len(binary_blocks) // (codewordlength)!= 0:
                    ValueError('The last codeword is not the right length')
            binarycodewords[-1][len(bitsofaddedcounter)] = str(decimal_to_binary(number_of_zero_padding)).zfill(howmanybitsin)
            uniformnameforoutput = binarycodewords

        if params.inner_error_correction == 'ltcode':
            splitltcode = makeltcode(uniformnameforoutput, params.percent_of_symbols, params.index_carry_length, params.ltcode_header,howmanybitsin)
            return splitltcode
        else:
            return uniformnameforoutput

    def decimal_to_binary_list(self, decimal, length):
        # Convert the decimal number to a binary string with the specified length
        binary_str = bin(decimal)[2:].zfill(length)
        # Convert the binary string to a list of individual numerals
        binary_list = [int(digit) for digit in binary_str]
        return binary_list  


    def encode_linear_binom(self, data):
        """
        Encodes the provided data using binomial encoding.
        """
        try:
            oligos = self.params.library.library

            leftmotif = self.params.library.leftmotives
            rightmotif = self.params.library.rightmotives
            
            leftmotif = list(sorted(set(leftmotif)))
            rightmotif = list(sorted(set(rightmotif)))

            positionallibs = create_positional_libraries(self.params.library, len(leftmotif) + len(rightmotif) - 1)
            
            # create binary codewords
            binary_codewords = self.create_binary_codewords(data.data, self.params)

            # create combinatorial codes
            dnacombinatorialcodes = []

            # TODO: paralize this loop
            # translate the binary codewords to combinations of library oligos
            for codeword in range(len(binary_codewords)):
                possible_oligos = []
                for i in range(len(binary_codewords[codeword])):
                    
                    positionallibs[i + self.params.dna_barcode_length].sort()
                    leftmotifofpos = []
                    rightmotifofpos = []
                    for j in range(len(positionallibs[i+self.params.dna_barcode_length])):
                        leftmotifofpos.append(positionallibs[i+self.params.dna_barcode_length][j][:len(oligos[0])//2])
                        rightmotifofpos.append(positionallibs[i+self.params.dna_barcode_length][j][len(oligos[0])//2:])
                    leftmotifofpos.sort()
                    rightmotifofpos.sort()
                    leftmotifofpos = self.remove_complements(leftmotifofpos)
                    rightmotifofpos = self.remove_complements(rightmotifofpos)
                    leftmotifofpos = list(set(leftmotifofpos))
                    rightmotifofpos = list(set(rightmotifofpos))
                    leftmotifofpos.sort()
                    rightmotifofpos.sort()                
                    newpositionallibs =[''.join(combo) for combo in it.product(leftmotifofpos, rightmotifofpos)]


                    listofzeros = [0] * len(newpositionallibs)
                    decimal= int(binary_codewords[codeword][i], 2)
                    for j in range(self.params.sigma_amount):
                        listofzeros[j] = 1
                    counterofwhereblock = 0
                    looper = True
                    
                    while looper:
                        alluniqueoligos = []
                        dnaswithoutcomplimetns = []
                        for dnasinstep in range(len(listofzeros)):
                            if listofzeros[dnasinstep] == 1:
                                dnaswithoutcomplimetns.append(newpositionallibs[dnasinstep])
                        uniquemotivesleft = []
                        uniquemotivesright = []
                        for makinguniquecount in range(len(dnaswithoutcomplimetns)):
                            uniquemotivesleft.append(dnaswithoutcomplimetns[makinguniquecount][:len(oligos[0]) // 2])
                            uniquemotivesright.append(dnaswithoutcomplimetns[makinguniquecount][len(oligos[0]) // 2:])
                        
                        uniquemotivesleft = list(set(uniquemotivesleft))
                        uniquemotivesright = list(set(uniquemotivesright))
                        uniquecountleft = len(uniquemotivesleft)
                        uniquecountright = len(uniquemotivesright)
                        uniquemotivesleft.sort()
                        uniquemotivesright.sort()
                        for k in range(len(uniquemotivesleft)):
                            alluniqueoligos.append(uniquemotivesleft[k])
                        for k in range(len(uniquemotivesright)):
                            alluniqueoligos.append(uniquemotivesright[k])
                        
                        uniquecountoligos = uniquecountleft + uniquecountright
                        lastblockborder = counterofwhereblock
                        counterofwhereblock = counterofwhereblock + 2 ** uniquecountoligos
                        if decimal<counterofwhereblock:
                            looper = False
                        else: 
                            for chekforonetomove in range(len(listofzeros)-1, -1, -1):
                                if listofzeros[chekforonetomove] == 0 and listofzeros[chekforonetomove - 1] == 1:
                                    l = chekforonetomove -1
                                    listofzeros[l] = 0
                                    listofzeros[l + 1] = 1
                                    ones_to_move = listofzeros[l+2:]
                                    listofzeros[l+2:] = [0] * len(ones_to_move)
                                    listofzeros[l+2:] = [1] * ones_to_move.count(1) + [0] * ones_to_move.count(0)
                                    break

                    howdeepintotheblock = decimal - lastblockborder
                    whichmotivesarecompliments = self.decimal_to_binary_list(howdeepintotheblock, uniquecountoligos)
                    for turningcompliments in range(len(alluniqueoligos)):
                        if turningcompliments < uniquecountleft and whichmotivesarecompliments[turningcompliments] == 1:
                            for checkingcompliments in range(len(dnaswithoutcomplimetns)):
                                if dnaswithoutcomplimetns[checkingcompliments][:len(oligos[0])//2] == alluniqueoligos[turningcompliments]:
                                    dnaswithoutcomplimetns[checkingcompliments] = complementmap(dnaswithoutcomplimetns[checkingcompliments][:len(oligos[0])//2])+dnaswithoutcomplimetns[checkingcompliments][len(oligos[0])//2:]
                        if turningcompliments >= uniquecountleft and whichmotivesarecompliments[turningcompliments] == 1:
                            for checkingcompliments in range(len(dnaswithoutcomplimetns)):
                                if dnaswithoutcomplimetns[checkingcompliments][len(oligos[0])//2:] == alluniqueoligos[turningcompliments]:
                                    dnaswithoutcomplimetns[checkingcompliments] = dnaswithoutcomplimetns[checkingcompliments][:len(oligos[0])//2]+complementmap(dnaswithoutcomplimetns[checkingcompliments][len(oligos[0])//2:])

                    possible_oligos.append(dnaswithoutcomplimetns)    
                dnacombinatorialcodes.append(possible_oligos)
            
            # translate the combinatorial codewords to a nested list of DNA oligotrios
            oligotrios = []
            
            for i in range(len(dnacombinatorialcodes)):
                oligotriostemp0 = []
                counterdna = create_counter_list(self.params.dna_barcode_length, len(positionallibs[0]) - 1, i)
                for j in range(len(dnacombinatorialcodes[i]) + self.params.dna_barcode_length):
                    oligotriostemp1 = []
                    if j < self.params.dna_barcode_length:
                        
                        for k in range(len(dnacombinatorialcodes[i][j])):
                            oligotriostemp2 = []
                            if j == 0:
                                oligotriostemp2.append(leftmotif[j//2]+complementmap(positionallibs[j][counterdna[j]][len(oligos[0])//2:]))
                                oligotriostemp2.append(positionallibs[j][counterdna[j]])
                                oligotriostemp2.append(complementmap(positionallibs[j][counterdna[j]][:len(oligos[0])//2])+rightmotif[j//2])
                            else:
                                if j % 2 == 0:
                                    oligotriostemp2.append(leftmotif[j//2]+complementmap(positionallibs[j][counterdna[j]][len(oligos[0])//2:]))
                                    oligotriostemp2.append(positionallibs[j][counterdna[j]])
                                    oligotriostemp2.append(complementmap(positionallibs[j][counterdna[j]][:len(oligos[0])//2])+rightmotif[(j)//2])
                                else:
                                    oligotriostemp2.append(complementmap(positionallibs[j][counterdna[j]][:len(oligos[0])//2])+complementmap(rightmotif[(j+1)//2]))
                                    oligotriostemp2.append(positionallibs[j][counterdna[j]])
                                    oligotriostemp2.append(complementmap(leftmotif[(j-1)//2])+complementmap(positionallibs[j][counterdna[j]][len(oligos[0])//2:]))
                            
                            oligotriostemp1.append(oligotriostemp2)
                        
                    else:
                        for k in range(len(dnacombinatorialcodes[i][j-self.params.dna_barcode_length])):
                            oligotriostemp2 = []
                            if j == 0:
                                oligotriostemp2.append(leftmotif[j//2]+complementmap(dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k][len(oligos[0])//2:]))
                                oligotriostemp2.append(dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k])
                                oligotriostemp2.append(complementmap(dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k][:len(oligos[0])//2])+rightmotif[j//2])
                            else:
                                if j % 2 == 0:
                                    oligotriostemp2.append(leftmotif[j//2]+complementmap(dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k][len(oligos[0])//2:]))
                                    oligotriostemp2.append(dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k])
                                    oligotriostemp2.append(complementmap(dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k][:len(oligos[0])//2])+rightmotif[(j)//2])
                                else:
                                    oligotriostemp2.append(complementmap(dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k][:len(oligos[0])//2])+complementmap(rightmotif[(j+1)//2]))
                                    oligotriostemp2.append((dnacombinatorialcodes[i][j-self.params.dna_barcode_length][k]))
                                    oligotriostemp2.append(complementmap(leftmotif[(j - 1) // 2]) + complementmap(dnacombinatorialcodes[i][j - self.params.dna_barcode_length][k][len(oligos[0]) // 2:]))
                            
                            oligotriostemp1.append(oligotriostemp2)

                    oligotriostemp0.append(oligotriostemp1)
                oligotrios.append(oligotriostemp0)


            info = {
                "number_of_codewords": len(oligotrios),
                "length_of_each_codeword": [len(codeword) for codeword in oligotrios],
                "barcode_length": self.params.dna_barcode_length,
                "data_length": len(data.data)
            }

            if self.params.theory == 'yes':
                theory_output = []
                for i in range(len(oligotrios)):
                    theory_output.append(self.combine_elements_randomized(oligotrios[i]))
                return theory_output, info

        except Exception as e:
            if self.logger:
                self.logger.error(f"Error encoding data: {str(e)}")
                self.logger.error(traceback.format_exc())
            oligotrios = None
            info = None
        
        return oligotrios, info
        
