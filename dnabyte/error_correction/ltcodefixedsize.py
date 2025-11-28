import numpy as np
import math as m
import random
#from mi_dna_disc.logging_config import logger

def translate_and_join(numbers, length):
    binary_strings = [format(num, f'0{length}b') for num in numbers]
    joined_string = ''.join(binary_strings)
    return joined_string

def split_encoded_symbol(encoded_symbol, ninputlength, indexcarrylength):
    try:
        if len(encoded_symbol) < ninputlength + indexcarrylength:
            raise ValueError("The length of the binary string is less than the sum of ninputlength and indexcarrylength.")
    except:
        #logger.info('Decoding failed: The length of the binary string is less than the sum of ninputlength and indexcarrylength.', exc_info=True)
        pass

    howmanyarethere = encoded_symbol[:ninputlength]
    indices_part = encoded_symbol[-indexcarrylength:]
    main_encoded_symbol = encoded_symbol[ninputlength:-indexcarrylength] if ninputlength + indexcarrylength < len(encoded_symbol) else ''
    
    return howmanyarethere, main_encoded_symbol, indices_part

def split_binary_string_from_back(binary_string, n):

    return [binary_string[max(0, len(binary_string) - (i + 1) * n):len(binary_string) - i * n] for i in range((len(binary_string) + n - 1) // n)]

def remove_zeros_after_last_non_zero(lst):
    last_non_zero_index = -1
    for i in range(len(lst)):
        if lst[i] != 0:
            last_non_zero_index = i
    if last_non_zero_index == -1:
        return [0]
    return lst[:last_non_zero_index + 1]


def robust_soliton_distribution(N, c=0.1, delta=0.5, max_degree=None):
    if max_degree is None or max_degree > N:
        max_degree = N

    R = c * np.log(N / delta) * np.sqrt(N)
    p = [1 / N]
    for i in range(2, max_degree + 1):
        p.append(1 / (i * (i - 1)))
    p = np.array(p)
    p /= p.sum()

    # Add the extra term for the robust soliton distribution
    denominator = (np.arange(1, max_degree + 1) * (np.arange(1, max_degree + 1) - 1))
    denominator[0] = N
    tau = R/denominator
   

    # tau = R / (np.arange(1, max_degree + 1) * (np.arange(1, max_degree + 1) - 1))
    # tau[0] = R / N  # Adjust the first element to avoid division by zero
    p += tau
    p /= p.sum()
    
    # Adjust the distribution to ensure it sums to 1
    if max_degree < N:
        p = np.append(p, [0] * (N - max_degree))
        p /= p.sum()
    
    return p

def encode_lt(input_strings, num_symbols, indexcarrylength, ninputlength):
    ninput = len(input_strings)

    try:
        if ninput > 2 ** ninputlength:
            raise ValueError("The number of input strings is too large for the given input length.")
        else:
            ninputrember = bin(ninput)[2:].zfill(ninputlength)
    except:
        #logger.info('Encoding failed: The number of input strings is too large for the given input length.', exc_info=True)
        pass

    howmanybitsforoneindex = max(1,int(m.ceil(m.log2(ninput))))
    max_degree = m.floor(indexcarrylength / howmanybitsforoneindex)
    p = robust_soliton_distribution(ninput, max_degree=max_degree)
    encoded_symbols = []
    for _ in range(num_symbols):
        degree = np.random.choice(np.arange(1, ninput + 1), p=p)
        chosen_indices = random.sample(range(ninput), min(degree, max_degree))
        chosen_indices.sort(reverse=True)
        chosen_indices_binary = translate_and_join(chosen_indices, howmanybitsforoneindex)
        if len(chosen_indices_binary) < indexcarrylength:
            chosen_indices_binary = '0' * (indexcarrylength - len(chosen_indices_binary)) + chosen_indices_binary
        encoded_symbol = input_strings[chosen_indices[0]]
    
        for idx in chosen_indices[1:]:
            encoded_symbol = ''.join(str(int(a) ^ int(b)) for a, b in zip(encoded_symbol, input_strings[idx]))

        encoded_symbols.append(ninputrember + encoded_symbol + chosen_indices_binary)
        
    return encoded_symbols

def decode_lt(encoded_symbols, indexcarrylength, ninputlength):
    
    howmanyarethereonavragelist = []
    for symbols in encoded_symbols:
        howmanyarethere, main_encoded_symbol, indices_part = split_encoded_symbol(symbols, ninputlength, indexcarrylength)
        howmanyarethereonavragelist.append(howmanyarethere)

    if howmanyarethereonavragelist == []:
        raise ValueError("No encoded symbols to decode.")
    
    howmanyarethereonavrage = max(set(howmanyarethereonavragelist), key=howmanyarethereonavragelist.count)
    howmanyaretheredec = int(howmanyarethereonavrage, 2)
    decoded_strings = [None] * howmanyaretheredec
    unresolved_symbols = set(range(howmanyaretheredec))
    howmanybitsforoneindex = int(m.ceil(m.log2(howmanyaretheredec)))
    if howmanybitsforoneindex == 0:
        raise ValueError("LTcode header has been corrupted and cannot determine the original amount of messages.")
    
    encoding_struc = []
    for symbols in encoded_symbols:
        howmanyaretherenotuse, main_encoded_symbol, indices_part = split_encoded_symbol(symbols, ninputlength, indexcarrylength)
        indecesinlist = split_binary_string_from_back(indices_part, howmanybitsforoneindex)
        chosen_indicesa = remove_zeros_after_last_non_zero([int(idx, 2) for idx in indecesinlist])
        if all(chosen_indicesa[idx] < howmanyaretheredec for idx in range(len(chosen_indicesa))):
            
            encoding_struc.append((chosen_indicesa, main_encoded_symbol))

    try:
        while unresolved_symbols:
            progress = False
            new_encoded_symbols = []

            for chosen_indices, encoded_symbol in encoding_struc:
                if len(chosen_indices) == 1 and chosen_indices[0] in unresolved_symbols:
                    decoded_strings[chosen_indices[0]] = encoded_symbol
                    unresolved_symbols.remove(chosen_indices[0])
                    progress = True
                else:
                    for idx in chosen_indices:
                        if len(decoded_strings) < idx:
                            raise ValueError("Data to currupted to decode.")
                        if decoded_strings[idx] is not None:
                            encoded_symbol = ''.join(str(int(a) ^ int(b)) for a, b in zip(encoded_symbol, decoded_strings[idx]))
                    chosen_indices = [idx for idx in chosen_indices if decoded_strings[idx] is None]
                    if chosen_indices:
                        new_encoded_symbols.append((chosen_indices, encoded_symbol))

            encoding_struc = new_encoded_symbols

            if not progress:
                break
            
    except:
        #logger.info('Decoding failed: Too many errors in header of fountaincode in every codeword.', exc_info=True)
        return [x for x in decoded_strings if x is not None], False
    if None in decoded_strings:
        #logger.info('Decoding failed: Either too many errors in header of fountaincode in every codeword or too many lost codewords.')
        return [x for x in decoded_strings if x is not None], False

    return decoded_strings, True


