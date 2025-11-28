import random

class Oligo:
    """
        This class models a single or double stranded oligonucleotide. It is always instanciated as single stranded,
        from which double stranded oligonucleotides can be created by hybridising it with another oligo.

        Attributes:
            motifs (tuple): A pair of motifs where the first motif is lowercase (representing the 5' end) and the 
            second motif is uppercase (representing the 3' end). This attribute is None if no motifs are provided.
            type (str): Indicates the type of the oligonucleotide. Defaults to 'single_stranded'.

        Raises:
            ValueError: If `motifs` is provided but is not a valid motif pair (i.e., a tuple of two elements 
            where the first element is lowercase and the second is uppercase).

        Parameters:
            motifs (tuple, optional): A tuple containing a pair of motifs to initialize the oligonucleotide. 
            Defaults to None.
        """

    def __init__(self, motifs=None, sequence=None):
        self.motifs = motifs
        self.sequence = sequence

        if motifs:
            if isinstance(motifs[0], tuple):
                self.type = 'double_stranded'
            if isinstance(motifs[0], str):
                self.type = 'single_stranded'
        # TODO: This causes the oligo assignment not to work
        # Determine the type based on the structure of motifs or sequence
        # if motifs is not None:
        #     if isinstance(motifs, list):
        #         self.type = 'double_stranded'
        #     elif isinstance(motifs, tuple):
        #         self.type = 'single_stranded'
        #     else:
        #         raise ValueError("Invalid motif structure.")

        # if sequence is not None:
        #     if isinstance(sequence, list) and all(isinstance(s, str) for s in sequence):
        #         self.type = 'double_stranded'
        #     elif isinstance(sequence, str):
        #         self.type = 'single_stranded'
        #     else:
        #         raise ValueError("Invalid sequence structure.")
        # else:
        #     #raise ValueError("Oligo must be initialized with motifs or sequence.")
        #     return None

    def end(self, strand=None, side=None):
        # If only one argument is provided, it is considered as 'side'
        if strand is not None and side is None:
            side, strand = strand, None

        if self.type == 'single_stranded':
            if side == '5' or side == 'r':
                return self.motifs[1]
            elif side == '3' or side == 'l':
                return self.motifs[0]

        else:
            if strand == 'f' and (side == '5' or side == 'r'):
                return self.motifs[0][-1]
            elif strand == 'f' and (side == '3' or side == 'l'):
                return self.motifs[0][0]
            elif strand == 'r' and (side == '5' or side == 'l'):
                return self.motifs[1][0]
            elif strand == 'r' and (side == '3' or side == 'r'):
                return self.motifs[1][-1]


    def __str__(self):

        if self.motifs and hasattr(self, 'motifs') and self.type == 'single_stranded':
            formatted_motif_1 = f" {self.motifs[0]}  " if len(self.motifs[0]) == 1 else f" {self.motifs[0]} "
            formatted_motif_2 = f" {self.motifs[1]}  " if len(self.motifs[1]) == 1 else f" {self.motifs[1]} "
            return "3':  ( " + formatted_motif_1 + " - " + formatted_motif_2 + " )  :5' \n"

        elif self.motifs and hasattr(self, 'motifs') and self.type == 'double_stranded':
            output = ""
            for index, strand in enumerate(self.motifs):
                line = "3':  (" if index == 0 else "5':  ("
                
                for i, motif in enumerate(strand):
                    if motif:
                        line += f" {motif}  " if len(motif) == 1 else f" {motif} "
                    else:
                        line += "    "
                    if i < len(strand) - 1:  # Check if it's not the last motif
                        line += "-"
                line += ")  :5'  \n" if index == 0 else ")  :3' \n"
                output += line
            output += "\n"
            return output

        elif self.sequence and hasattr(self, 'sequence') and self.type == 'single_stranded':
            return "3': " + self.sequence + " :5'"

        elif self.sequence and hasattr(self, 'sequence') and self.type == 'double_stranded':
            return "3': " + self.sequence[0] + " :5' \n5': " + self.sequence[1] + " :3'"

        else:
            return "Oligo does not have properly defined motifs."
        
    wc_pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


    def mirror(self):
        """
        Returns the identical oligonucleotide, where forward and reverse strand are exchanged.
        """
        if self.type == 'single_stranded':
            return self
        else: 
            return Oligo(motifs=(self.motifs[1][::-1], self.motifs[0][::-1]))




def create_positional_libraries(library, codewordlength):
    positionallibaraies = []
    leftmotif = library.leftmotives
    rightmotif = library.rightmotives
    leftmotif = list(set(leftmotif))
    rightmotif = list(set(rightmotif))
    leftmotif.sort()
    rightmotif.sort()
    oligos = library.library
    oligos.sort()

    for i in range(codewordlength):
        positionallibrary = []
        if i == 0:
            positionallibrary=[x for x in oligos if leftmotif[i//2] not in x and complementmap(leftmotif[i//2]) not in x and rightmotif[i//2] not in x and complementmap(rightmotif[i//2]) not in x]

        else:
            if i%2==0:
                left_motif = leftmotif[(i)//2]
                right_motif = rightmotif[(i)//2]
                left_motif_compliment = complementmap(left_motif)
                right_motif_compliment = complementmap(right_motif)

                positionallibrary = [
                    x for x in oligos
                    if (left_motif not in x and left_motif_compliment not in x and
                        right_motif not in x and right_motif_compliment not in x)]
            else:                
                left_motif = leftmotif[(i-1)//2]
                right_motif = rightmotif[(i+1)//2]
                left_motif_compliment = complementmap(left_motif)
                right_motif_compliment = complementmap(right_motif)

                positionallibrary = [
                    x for x in oligos
                    if (left_motif not in x and left_motif_compliment not in x and
                        right_motif not in x and right_motif_compliment not in x)
                ]
        positionallibrary.sort()
        positionallibaraies.append(positionallibrary)
    return positionallibaraies
   

def complement(motif, motif_pair):
    # Define Watson-Crick pairing
    
    #this dictionary should be loaded from the library.
    # motif_pair = {'a': 'a*', 'b': 'b*', 'c': 'c*', 'd': 'd*', 'e': 'e*', 'f': 'f*', 'g': 'g*', 'h': 'h*',
    #               'a*': 'a', 'b*': 'b', 'c*': 'c', 'd*': 'd', 'e*': 'e', 'f*': 'f', 'g*': 'g', 'h*': 'h',
    #               'A': 'A*', 'B': 'B*', 'C': 'C*', 'D': 'D*', 'E': 'E*', 'F': 'F*', 'G': 'G*', 'H': 'H*',
    #               'A*': 'A', 'B*': 'B', 'C*': 'C', 'D*': 'D', 'E*': 'E', 'F*': 'F', 'G*': 'G', 'H*': 'H', 
    #               'empty':'empty*', None: None, 'barcode':'barcode*', 'barcode*':'barcode',}
        
    if (isinstance(motif, str) or motif == None) and motif in motif_pair:
        return motif_pair[motif]
    
    elif (isinstance(motif, tuple) or motif == None) and all(m in motif_pair for m in motif):
        return tuple([motif_pair[m] for m in motif])
    
    else:
        raise ValueError("Invalid motif.")
    

def nucleotide_complement(nucleotide):
    # Define Watson-Crick pairing
    nucleotide_pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    if isinstance(nucleotide, str):
        if nucleotide in nucleotide_pairs:
            return nucleotide_pairs[nucleotide]
        
    elif isinstance(nucleotide, list):
        return [nucleotide_pairs[nuc] for nuc in nucleotide]

    else:
        raise ValueError("Invalid data type.")      
    

def complementmap(string):
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


def translate_element(element, translation_dictleft,translation_dictright):
    
    if isinstance(element, list):
        return [translate_element(sub_element,translation_dictleft,translation_dictright) for sub_element in element]
    else:
        return Oligo([translation_dictleft[element[:len(element)//2]],translation_dictright[element[len(element)//2:]]])
    

def translate_nested_list(nested_list, translation_dictleft,translation_dictright):
    return translate_element(nested_list, translation_dictleft,translation_dictright)

def process_tuple_list(tuple_list):
    def flip_tuple(t):
        return t[::-1]

    # Ensure the tuple containing 'A' is the first tuple
    if 'A' == tuple_list[1][0] or 'A' == tuple_list[1][-1]:
        tuple_list = (tuple_list[1], tuple_list[0])

    # Ensure 'A' is the first entry of the first tuple
    if tuple_list[0][0] != 'A':
        tuple_list = (flip_tuple(tuple_list[0]), flip_tuple(tuple_list[1]))

    return tuple_list

def process_tuple_list_empy(tuple_list):
    def flip_tuple(t):
        return t[::-1]

    # Ensure the tuple containing 'A' is the first tuple
    if 'empty' == tuple_list[1][0] or 'empty' == tuple_list[1][-1]:
        tuple_list = (tuple_list[1], tuple_list[0])

    # Ensure 'A' is the first entry of the first tuple
    if tuple_list[0][0] != 'empty':
        tuple_list = (flip_tuple(tuple_list[0]), flip_tuple(tuple_list[1]))

    return tuple_list

def process_tuple_list_poly(tuple_list):
    def flip_tuple(t):
        return t[::-1]

    # Ensure the tuple containing 'A' is the first tuple
    if 'c0' == tuple_list[1][0] or 'c0' == tuple_list[1][-1]:
        tuple_list = (tuple_list[1], tuple_list[0])

    # Ensure 'A' is the first entry of the first tuple
    if tuple_list[0][0] != 'c0':
        tuple_list = (flip_tuple(tuple_list[0]), flip_tuple(tuple_list[1]))

    return tuple_list

def back_translate_nested_list(element, translation_dictleft,translation_dictright):
    # TODO. fix header
    retranslated = []
    translation_dictleft = {v: k for k, v in translation_dictleft.items()}
    translation_dictright = {v: k for k, v in translation_dictright.items()}
    


    for i in range(len(element)):
        perpool = []
        for j in range(len(element[i])):
            
            fiveprime = []
            threeprime = []
            correcterdlist = process_tuple_list(element[i][j].pool[0].motifs)
            if type(correcterdlist[0]) == str or type(correcterdlist[1]) == str:
                pass
            else:
                for k in range(len(correcterdlist[0])):
                    if correcterdlist[0][k] in translation_dictleft.keys() or correcterdlist[1][k] in translation_dictleft.keys():
                        if correcterdlist[0][k] != None:
                            fiveprime.append(translation_dictleft[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            threeprime.append(translation_dictleft[correcterdlist[1][k]])
                    else:
                        if correcterdlist[0][k] != None:
                            fiveprime.append(translation_dictright[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            threeprime.append(translation_dictright[correcterdlist[1][k]])
                
                perpool.append([''.join(fiveprime),''.join(threeprime)])
        retranslated.append(perpool)
    return retranslated

def back_translate_nested_list_real(element, translation_dictleft,translation_dictright):
    # TODO. fix header
    retranslated = []
    translation_dictleft = {v: k for k, v in translation_dictleft.items()}
    translation_dictright = {v: k for k, v in translation_dictright.items()}

    
    # if len(element[0][0].pool) != 1: 
    #         ValueError("This did not perfectly hybridised")
    perpool = []
    for i in range(len(element.pool)):
        if 'A' == element.pool[i].motifs[1][0] or 'A' == element.pool[i].motifs[1][-1] or 'A' == element.pool[i].motifs[0][0] or 'A' == element.pool[i].motifs[0][-1]:
            
            fiveprime = []
            threeprime = []
            correcterdlist = process_tuple_list(element.pool[i].motifs)
            if type(correcterdlist[0]) == str or type(correcterdlist[1]) == str:
                pass
            else:
                for k in range(len(correcterdlist[0])):
                    if correcterdlist[0][k] in translation_dictleft.keys() or correcterdlist[1][k] in translation_dictleft.keys():
                        if correcterdlist[0][k] != None:
                            fiveprime.append(translation_dictleft[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            threeprime.append(translation_dictleft[correcterdlist[1][k]])
                    else:
                        if correcterdlist[0][k] != None:
                            fiveprime.append(translation_dictright[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            threeprime.append(translation_dictright[correcterdlist[1][k]])
                
                perpool.append([''.join(fiveprime),''.join(threeprime)])
        else:
            pass
    return perpool

def back_translate_nested_list_chain_real(element, translation_dictleft,translation_dictright):
    # TODO. fix header
    retranslated = []
    translation_dictleft = {v: k for k, v in translation_dictleft.items()}
    translation_dictright = {v: k for k, v in translation_dictright.items()}
    perpool = []
    for i in range(len(element.pool)):
        
        fiveprime = []
        threeprime = []
        if 'empty' == element.pool[i].motifs[1][0] or 'empty' == element.pool[i].motifs[1][-1] or 'empty' == element.pool[i].motifs[0][0] or 'empty' == element.pool[i].motifs[0][-1]:
            correcterdlist = process_tuple_list_empy(element.pool[i].motifs)

            if type(correcterdlist[0]) == str or type(correcterdlist[1]) == str:
                pass   
            else:         
                for k in range(len(correcterdlist[0])):
                    if correcterdlist[0][k] in translation_dictleft.keys() or correcterdlist[1][k] in translation_dictleft.keys():
                        if correcterdlist[0][k] != None:
                            if correcterdlist[0][k] != 'empty':
                                fiveprime.append(translation_dictleft[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            if correcterdlist[1][k] != 'empty':
                                threeprime.append(translation_dictleft[correcterdlist[1][k]])
                    else:
                        if correcterdlist[0][k] != None:
                            if correcterdlist[0][k] != 'empty':
                                fiveprime.append(translation_dictright[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            if correcterdlist[1][k] != 'empty':
                                threeprime.append(translation_dictright[correcterdlist[1][k]])
                    
                perpool.append([''.join(fiveprime),''.join(threeprime)])
        else:
            pass
    return perpool

def back_translate_nested_list_chain(element, translation_dictleft,translation_dictright):
    # TODO. fix header
    retranslated = []
    translation_dictleft = {v: k for k, v in translation_dictleft.items()}
    translation_dictright = {v: k for k, v in translation_dictright.items()}
    

    perpool = []
    for i in range(len(element)):
        fiveprime = []
        threeprime = []
        correcterdlist = process_tuple_list_empy(element[i].pool[0].motifs)

        if type(correcterdlist[0]) == str or type(correcterdlist[1]) == str:
            pass
        else:
            for k in range(len(correcterdlist[0])):
                if correcterdlist[0][k] in translation_dictleft.keys() or correcterdlist[1][k] in translation_dictleft.keys():
                    if correcterdlist[0][k] != None:
                        if correcterdlist[0][k] != 'empty':
                            fiveprime.append(translation_dictleft[correcterdlist[0][k]])
                    if correcterdlist[1][k] != None:
                        if correcterdlist[1][k] != 'empty':
                            threeprime.append(translation_dictleft[correcterdlist[1][k]])
                else:
                    if correcterdlist[0][k] != None:
                        if correcterdlist[0][k] != 'empty':
                            fiveprime.append(translation_dictright[correcterdlist[0][k]])
                    if correcterdlist[1][k] != None:
                        if correcterdlist[1][k] != 'empty':
                            threeprime.append(translation_dictright[correcterdlist[1][k]])
                
            perpool.append([''.join(fiveprime),''.join(threeprime)])
    return perpool

def translate_element_poly(element, library):
    

    if isinstance(element, list):
        return [translate_element_poly(sub_element, library) for sub_element in element]
    else:
        # NOTE: there is an edgecase where this doesn't work
        if len(element) == len(library.messages[0]):
            if element[:len(element)//2] in library.messageleft:
                return Oligo([library.messageleft[element[:len(element)//2]],library.messageright[element[len(element)//2:]]])
            elif element[:len(element)//2] in library.messageright:
                return Oligo([library.messageright[element[:len(element)//2]],library.messageleft[element[len(element)//2:]]])
        elif len(element) == len(library.generic[0])+len(library.messages[0])//2:
            
            if element[:len(library.generic[0])] in library.generic or complementmap(element[:len(library.generic[0])]) in library.generic:
                if element[len(library.generic[0]):] in library.messageleft:
                    return Oligo([library.genericlib[element[:len(library.generic[0])]],library.messageleft[element[len(library.generic[0]):]]])
                else:
                    return Oligo([library.genericlib[element[:len(library.generic[0])]],library.messageright[element[len(library.generic[0]):]]])
                
            elif element[len(library.messages[0])//2:] in library.generic or complementmap(element[len(library.messages[0])//2:]) in library.generic:
                if element[:len(library.messages[0])//2] in library.messageleft:
                    return Oligo([library.messageleft[element[:len(library.messages[0])//2]],library.genericlib[element[len(library.messages[0])//2:]]])
                else:
                    return Oligo([library.messageright[element[:len(library.messages[0])//2]],library.genericlib[element[len(library.messages[0])//2:]]])
               
        elif len(element) == len(library.generic[0])+len(library.position[0]):
            if element[:len(library.generic[0])] in library.generic or complementmap(element[:len(library.generic[0])]) in library.generic:
                return Oligo([library.genericlib[element[:len(library.generic[0])]],library.connectorlib[element[len(library.generic[0]):]]])
                
            elif element[len(library.position[0]):] in library.generic or complementmap(element[len(library.position[0]):]) in library.generic:
                return Oligo([library.connectorlib[element[:len(library.position[0])]],library.genericlib[element[len(library.position[0]):]]])
           
        else:
            ValueError('something went wrong')

        

def translate_nested_list_poly(nested_list, libraray):
    return translate_element_poly(nested_list, libraray)


def back_translate_nested_list_poly(element, library):
    genericlib = library.genericlib
    genericlib = {v: k for k, v in genericlib.items()}
    connectorlib = library.connectorlib
    connectorlib = {v: k for k, v in connectorlib.items()}
    messageleft = library.messageleft
    messageleft = {v: k for k, v in messageleft.items()}
    messageright = library.messageright
    messageright = {v: k for k, v in messageright.items()}

    perpool = []
    for i in range(len(element)):
        fiveprime = []
        threeprime = []

        
        
        correcterdlist = process_tuple_list_poly(element[i].pool[0].motifs)
        if type(correcterdlist[0]) == str or type(correcterdlist[1]) == str:
            pass
        else:
            for k in range(len(correcterdlist[0])):
                if correcterdlist[0][k] in genericlib.keys() or correcterdlist[1][k] in genericlib.keys():
                    if correcterdlist[0][k] != None:
                        if correcterdlist[0][k] != 'empty':
                            fiveprime.append(genericlib[correcterdlist[0][k]])
                    if correcterdlist[1][k] != None:
                        if correcterdlist[1][k] != 'empty':
                            threeprime.append(genericlib[correcterdlist[1][k]])
                elif correcterdlist[0][k] in connectorlib.keys() or correcterdlist[1][k] in connectorlib.keys():
                    if correcterdlist[0][k] != None:
                        if correcterdlist[0][k] != 'empty':
                            fiveprime.append(connectorlib[correcterdlist[0][k]])
                    if correcterdlist[1][k] != None:
                        if correcterdlist[1][k] != 'empty':
                            threeprime.append(connectorlib[correcterdlist[1][k]])
                elif correcterdlist[0][k] in messageleft.keys() or correcterdlist[1][k] in messageleft.keys():
                    if correcterdlist[0][k] != None:
                        if correcterdlist[0][k] != 'empty':
                            fiveprime.append(messageleft[correcterdlist[0][k]])
                    if correcterdlist[1][k] != None:
                        if correcterdlist[1][k] != 'empty':
                            threeprime.append(messageleft[correcterdlist[1][k]])
                else:
                    if correcterdlist[0][k] != None:
                        if correcterdlist[0][k] != 'empty':
                            fiveprime.append(messageright[correcterdlist[0][k]])
                    if correcterdlist[1][k] != None:
                        if correcterdlist[1][k] != 'empty':
                            threeprime.append(messageright[correcterdlist[1][k]])
                
            perpool.append([''.join(fiveprime),''.join(threeprime)])
    return perpool

def back_translate_nested_list_poly_binom(element, library):
    backtranslated = []
    for i in range(len(element)):
        backtranslated.append(back_translate_nested_list_poly(element[i], library)	)

    return backtranslated
    
def back_translate_nested_list_poly_binom_real(element, library):
    genericlib = library.genericlib
    genericlib = {v: k for k, v in genericlib.items()}
    connectorlib = library.connectorlib
    connectorlib = {v: k for k, v in connectorlib.items()}
    messageleft = library.messageleft
    messageleft = {v: k for k, v in messageleft.items()}
    messageright = library.messageright
    messageright = {v: k for k, v in messageright.items()}

    perpool = []
    
    for i in range(len(element.pool)):
        fiveprime = []
        threeprime = []
        
        if 'c0' == element.pool[i].motifs[1][0] or 'c0' == element.pool[i].motifs[1][-1] or 'c0' == element.pool[i].motifs[0][0] or 'c0' == element.pool[i].motifs[0][-1]:
            correcterdlist = process_tuple_list_poly(element.pool[i].motifs)
            if type(correcterdlist[0]) == str or type(correcterdlist[1]) == str:
                pass
            else:
                for k in range(len(correcterdlist[0])):
                    if correcterdlist[0][k] in genericlib.keys() or correcterdlist[1][k] in genericlib.keys():
                        if correcterdlist[0][k] != None:
                            if correcterdlist[0][k] != 'empty':
                                fiveprime.append(genericlib[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            if correcterdlist[1][k] != 'empty':
                                threeprime.append(genericlib[correcterdlist[1][k]])
                    elif correcterdlist[0][k] in connectorlib.keys() or correcterdlist[1][k] in connectorlib.keys():
                        if correcterdlist[0][k] != None:
                            if correcterdlist[0][k] != 'empty':
                                fiveprime.append(connectorlib[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            if correcterdlist[1][k] != 'empty':
                                threeprime.append(connectorlib[correcterdlist[1][k]])
                    elif correcterdlist[0][k] in messageleft.keys() or correcterdlist[1][k] in messageleft.keys():
                        if correcterdlist[0][k] != None:
                            if correcterdlist[0][k] != 'empty':
                                fiveprime.append(messageleft[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            if correcterdlist[1][k] != 'empty':
                                threeprime.append(messageleft[correcterdlist[1][k]])
                    else:
                        if correcterdlist[0][k] != None:
                            if correcterdlist[0][k] != 'empty':
                                fiveprime.append(messageright[correcterdlist[0][k]])
                        if correcterdlist[1][k] != None:
                            if correcterdlist[1][k] != 'empty':
                                threeprime.append(messageright[correcterdlist[1][k]])
                    
                perpool.append([''.join(fiveprime),''.join(threeprime)])
        else:
            pass
        
    return perpool

