import math
import numpy as np
from random import choices
import random

class Symbol:
    def __init__(self, index, degree, data,neighbors):
        self.index = index
        self.degree = degree
        self.data = data
        self.neighbors=neighbors


def ideal_distribution(N):
    """ Create the ideal soliton distribution. 
    In practice, this distribution gives not the best results
    Cf. https://en.wikipedia.org/wiki/Soliton_distribution
    """
    EPSILON = 0.0001

    probabilities = [0, 1 / N]
    probabilities += [1 / (k * (k - 1)) for k in range(2, N+1)]
    probabilities_sum = sum(probabilities)

    assert probabilities_sum >= 1 - EPSILON and probabilities_sum <= 1 + EPSILON, "The ideal distribution should be standardized"
    return probabilities


def robust_distribution(N):
    """ Create the robust soliton distribution. 
    This fixes the problems of the ideal distribution
    Cf. https://en.wikipedia.org/wiki/Soliton_distribution
    """
    ROBUST_FAILURE_PROBABILITY = 0.001
    EPSILON = 0.0001
    # The choice of M is not a part of the distribution ; it may be improved
    # We take the median and add +1 to avoid possible division by zero 
    M = N // 2 + 1 
    R = N / M

    extra_proba = [0] + [1 / (i * M) for i in range(1, M)]
    extra_proba += [math.log(R / ROBUST_FAILURE_PROBABILITY) / M]  # Spike at M
    extra_proba += [0 for k in range(M+1, N+1)]

    probabilities = np.add(extra_proba, ideal_distribution(N))
    probabilities /= np.sum(probabilities)
    probabilities_sum = np.sum(probabilities)

    assert probabilities_sum >= 1 - EPSILON and probabilities_sum <= 1 + EPSILON, "The robust distribution should be standardized"
    return probabilities

def get_degrees_from(distribution_name, N, k):
    """ Returns the random degrees from a given distribution of probabilities.
    The degrees distribution must look like a Poisson distribution and the 
    degree of the first drop is 1 to ensure the start of decoding.
    """

    if distribution_name == "ideal":
        probabilities = ideal_distribution(N)
    elif distribution_name == "robust":
        probabilities = robust_distribution(N)
    else:
        probabilities = None
    
    population = list(range(0, N+1))
    return [1] + choices(population, probabilities, k=k-1)

def generate_indexes(symbol_index, degree, blocks_quantity):
    """Randomly get `degree` indexes, given the symbol index as a seed

    Generating with a seed allows saving only the seed (and the amount of degrees) 
    and not the whole array of indexes. That saves memory, but also bandwidth when paquets are sent.

    The random indexes need to be unique because the decoding process uses dictionnaries for performance enhancements.
    Additionnally, even if XORing one block with itself among with other is not a problem for the algorithm, 
    it is better to avoid uneffective operations like that.

    To be sure to get the same random indexes, we need to pass 
    """
    random.seed(symbol_index)
    indexes = random.sample(range(blocks_quantity), degree)

    return indexes, degree

def encode(blocks, drops_quantity):
    blocks_n = len(blocks)
    assert blocks_n <= drops_quantity, "Because of the unicity in the random neighbors, it is need to drop at least the same amount of blocks"
    # Generate random indexes associated to random degrees, seeded with the symbol id
    random_degrees = get_degrees_from("robust", blocks_n, k=drops_quantity)
    smbolsend=[]
    for i in range(drops_quantity):
        # Get the random selection, generated precedently (for performance)
        selection_indexes, deg = generate_indexes(i, random_degrees[i], blocks_n)
        drop = blocks[selection_indexes[0]]
        nigh=[]
        nigh.append(selection_indexes[0])
        for n in range(1, deg):
            drop = np.bitwise_xor(drop, blocks[selection_indexes[n]])
            nigh.append(selection_indexes[n])
            # drop = drop ^ blocks[selection_indexes[n]] # according to my tests, this has the same performance

        # Create symbol, then log the process
        dropbin=np.base_repr(int(drop))
        dropbin= '0' * (6 - len(dropbin)) + dropbin
        dropdna=translate_binary(dropbin)
        symbol = Symbol(index=i, degree=deg, data=dropdna,neighbors=nigh)
        smbolsend.append(symbol)

    return smbolsend

def translate_dna(dna_string):
    # Create a mapping from DNA bases to numbers
    dna_mapping = {'A': '0', 'C': '1', 'G': '0', 'T': '1'}

    # Translate the DNA string into numbers
    binary_string = ''.join(dna_mapping[dna] for dna in dna_string)

    # Represent the binary string as a binary number
    binary_number = int(binary_string, 2)
    base10 = np.base_repr(binary_number, base=10)
    return base10

def translate_binary(binary_string):
    # Create a mapping from binary digits to DNA bases
    binary_mapping = {0: ['A', 'G'], 1: ['C', 'T']}

    # Translate the binary string into DNA bases
    dna_string = ''.join(binary_mapping[int(digit)][i % 2] for i, digit in enumerate(binary_string))

    return dna_string

def dnafountaincode(data, drops_quantity):
    bindata=[]
    for i in data:
        bindata.append(int(translate_dna(i)))
    encoded_data = encode(bindata, drops_quantity)
    return encoded_data
