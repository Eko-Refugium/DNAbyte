from dnabyte.error_correction.error_correction_reed_solomon import ReedSolomon

def bitstream_to_bitmessages(bitstream, message_length, q):
    message_list = [[bitstream[i+j:i+j+q] for i in range(0, message_length*q, q)] for j in range(0, len(bitstream), message_length*q)]
    while '' in message_list[-1]:
        message_list[-1].remove('')
    return message_list

def bitmessages_to_bitstream(message_list):
    return ''.join([''.join(sublist) for sublist in message_list])

def bitmessages_to_decimalmessages(message_list):
    return [[int(bitstream, 2) for bitstream in sublist] for sublist in message_list]

def decimalmessages_to_bitmessages(message_list):
    return [[format(decimal, '08b') for decimal in sublist] for sublist in message_list]

def reed_solomon_encoding(message_list, p,q, message_length):
    coder=ReedSolomon(p,q, message_length)
    encoded_message_list = [coder.generate(message) for message in message_list]
    return encoded_message_list

def reed_solomon_decoding(encoded_message_list, p,q, message_length):
    coder=ReedSolomon(p,q, message_length)
    
    decoded_message_list=[]
    for i in range(len(encoded_message_list)):
        decoded_message_list.append(coder.correct(encoded_message_list[i]) )
    
    return decoded_message_list

def decimal_to_quaternary(decimal_number,q):
    if decimal_number == 0:
        return '0000'
    quaternary = ''
    for i in range(q//2):
        quaternary = str(decimal_number % 4) + quaternary
        decimal_number //= 4
    return quaternary

def quaternary_to_decimal(quaternary_number):
    decimal_number = 0
    for digit in quaternary_number:
        decimal_number = decimal_number * 4 + int(digit)
    return decimal_number

def decimalmessages_to_quartmessages(message_list,q):
    return [[decimal_to_quaternary(decimal,q) for decimal in sublist] for sublist in message_list]

def quartmessages_to_decimalmessages(message_list):
    return [[quaternary_to_decimal(tritmessage) for tritmessage in sublist] for sublist in message_list]

def join_list_of_lists(list_of_lists, separator=''):
    return separator.join([''.join(sublist) for sublist in list_of_lists])

def split_string_to_list_of_lists(quartmessage,quartmessage_length, bitlength):
    message_list = [[quartmessage[i+j:i+j+bitlength//2] for i in range(0, quartmessage_length*(bitlength//2), bitlength//2)] for j in range(0, len(quartmessage), quartmessage_length*(bitlength//2))]
    while '' in message_list[-1]:
        message_list[-1].remove('')
    #while [] in message_list[-1]:
     #   message_list[-1].remove('')
    return message_list

def quaternary_to_DNA(quaternarystring):
    mapping = {'0': 'A', '1': 'T', '2': 'C', '3': 'G'}
    return ''.join(mapping[digit] for digit in quaternarystring)

def DNA_to_quaternary(DNA_string):
    mapping = {'A': '0', 'T': '1', 'C': '2', 'G': '3'}
    return ''.join(mapping[letter] for letter in DNA_string)


def Err_decoding_Unconstrained(Dnaencoding,p,q, message_length):
    quartmessages = DNA_to_quaternary(Dnaencoding)
    splitlist = split_string_to_list_of_lists(quartmessages,p**q-1,q)
    decimalmessages = quartmessages_to_decimalmessages(splitlist)
    decoded_messages = []
    for i in range(len(decimalmessages)):
        lastbits=p**q-1-message_length
        final = []
        final = decimalmessages[i]
        for il in range(lastbits):
            final.pop()
        decoded_messages.append(final)
        
    bitmessages = decimalmessages_to_bitmessages(decoded_messages)
    bitstream = bitmessages_to_bitstream(bitmessages)
    return bitstream
