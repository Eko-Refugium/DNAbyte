from reedsolo import RSCodec
from dnabyte.error_correction.ltcodefixedsize import encode_lt, decode_lt

def bitstring_to_bytearray(bitstring):
    """
    Converts a bit string (e.g., '11001010') into a bytearray.

    Args:
        bitstring (str): A string of '0's and '1's representing bits.

    Returns:
        bytearray: The corresponding bytearray representation.
    """
    if len(bitstring) % 8 != 0:
        raise ValueError("Bit string length must be a multiple of 8")

    # Group the bit string into chunks of 8 bits and convert to bytes
    byte_array = bytearray(int(bitstring[i:i+8], 2) for i in range(0, len(bitstring), 8))
    return byte_array


def bytearray_to_bitstring(byte_array):
    """
    Converts a bytearray into a bit string (e.g., '11001010').

    Args:
        byte_array (bytearray): The bytearray to convert.

    Returns:
        str: A string of '0's and '1's representing the bits.
    """
    return ''.join(format(byte, '08b') for byte in byte_array)



def MakeReedSolomonCode(data,codewordlength,bitsinpos):
    
    reedsolomonlength = (codewordlength)*bitsinpos//8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    for codewords in data:
        bytearr = bitstring_to_bytearray(''.join(codewords))
        codeword = rs.encode(bytearr)
        reedsolobitstring = bytearray_to_bitstring(codeword)
        sploited = [reedsolobitstring[i:i+bitsinpos] for i in range(0, len(reedsolobitstring), bitsinpos)]
        reedsolomonencodedwords.append(sploited)
    
    return reedsolomonencodedwords

def MakeReedSolomonCodeSynthesis(data, errorlength):
    reedsolomonlength = (errorlength) // 8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    for codewords in data:
        bytearr = bitstring_to_bytearray(codewords)

        codewordrs = rs.encode(bytearr)
        reedsolobitstring = bytearray_to_bitstring(codewordrs)
        reedsolomonencodedwords.append(reedsolobitstring)
        
    return reedsolomonencodedwords

def undoreedsolomonsynthesis(data,errorlength):
    reedsolomonlength = (errorlength)//8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    valuechack = True
    for codewords in data:
        try:
            bytearr = bitstring_to_bytearray(codewords)
            decoded_bytearr = rs.decode(bytearr)
            decoded_bitstring = bytearray_to_bitstring(decoded_bytearr[0])
            reedsolomonencodedwords.append(decoded_bitstring)
        except:
            #logger.debug(f'Error in Reed Solomon decoding: codeword={codewords}, errorlength={errorlength}')
            valuechack = False
            reedsolomonencodedwords.append(codewords[:len(codewords)-errorlength])
    return reedsolomonencodedwords, valuechack

def undoreedsolomon(data,codewordlength,bitsinpos):
    
    reedsolomonlength = (codewordlength)*bitsinpos//8
    rs = RSCodec(reedsolomonlength)
    reedsolomonencodedwords = []
    valuechack = True
    for codewords in data:
        # try:
        binary_list = [bin(x)[2:].zfill(bitsinpos) for x in codewords]
        bytearr = bitstring_to_bytearray(''.join(binary_list))
        decoded_bytearr = rs.decode(bytearr)
        decoded_bitstring = bytearray_to_bitstring(decoded_bytearr[0])
        sploited = [decoded_bitstring[i:i+bitsinpos] for i in range(0, len(decoded_bitstring), bitsinpos)]
        decimal_list = [int(b, 2) for b in sploited]
        reedsolomonencodedwords.append(decimal_list)
        # except:
        #     #logger.debug(f'Error in Reed Solomon decoding: codeword={codewords}, errorlength={bitsinpos}')
        #     valuechack = False
        #     reedsolomonencodedwords.append(codewords[len(codewords)-bitsinpos:])
    
    return reedsolomonencodedwords, valuechack


def makeltcode(data, num_symbols, indexcarrylength, ninputlength,howmanybitsin):
    
    makingltcode = []
    for codewords in data:
        makingltcode.append(''.join(codewords))
    
    bitsofindexcarry = indexcarrylength*howmanybitsin
    encodedltcode = encode_lt(makingltcode, len(makingltcode)+int(len(makingltcode)*num_symbols),bitsofindexcarry,ninputlength*howmanybitsin )
    splitltcode = []
    for stringies in encodedltcode:
        splitltcode.append([stringies[i:i+howmanybitsin] for i in range(0, len(stringies), howmanybitsin)])
    return splitltcode

def makeltcodesynth(data, num_symbols, indexcarrylength, ninputlength, howmanybitsin):
    bitsofindexcarry = indexcarrylength * howmanybitsin
    encodedltcode = encode_lt(data, len(data) + int(len(data) * num_symbols), bitsofindexcarry, ninputlength * howmanybitsin)
    
    return encodedltcode

def undoltcode(data, indexcarrylength, ninputlength,howmanybitsin):
    encoded_symbols = []
    for codewords in data:
        binary_list = [bin(x)[2:].zfill(howmanybitsin) for x in codewords]
        binarystring = ''.join(binary_list)
        encoded_symbols.append(binarystring)
    bitsofindexcarry = indexcarrylength*howmanybitsin
    decodedltbin, checkforvalidlt=decode_lt(encoded_symbols,bitsofindexcarry,ninputlength*howmanybitsin)
    finischeddeclist = []
     
    for decodedlt in decodedltbin:
        decodedltsplit = [decodedlt[i:i+howmanybitsin] for i in range(0, len(decodedlt), howmanybitsin)]
        decimal_list = [int(b, 2) for b in decodedltsplit]
        finischeddeclist.append(decimal_list)
    return finischeddeclist, checkforvalidlt


def undoltcodesynth(data, indexcarrylength, ninputlength,howmanybitsin):
    
    bitsofindexcarry = indexcarrylength * howmanybitsin
    decodedltbin, checkforvalidlt = decode_lt(data, bitsofindexcarry, ninputlength * howmanybitsin)
    
    
    return decodedltbin, checkforvalidlt

