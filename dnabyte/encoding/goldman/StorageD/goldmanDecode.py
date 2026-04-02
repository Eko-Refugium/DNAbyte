# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 10:37:21 2021

@author: Wei QWiang
"""

import argparse
from tqdm import tqdm

def linuxCommand():
    parser = argparse.ArgumentParser()
    parser.add_argument("i", help = "input path")
    parser.add_argument("o", help = "output path")
    parser.add_argument("il", help = "index length", type=int)
    parser.add_argument("al", help = "add length", type = int)
    args = parser.parse_args()
    
    return args.i, args.o, args.il, args.al


def decodeNt(nt_seq):
    rotate_codes_dic = {'A': ['C', 'G', 'T'], 'C': ['G', 'T', 'A'], 'G': ['T', 'A', 'C'], 'T': ['A', 'C', 'G']}
    last_nt = "A"

    huffman_str = ""
    for nt in nt_seq:
        choose_list = rotate_codes_dic[last_nt]
        huffman_str += str(choose_list.index(nt))
        last_nt = nt

    return huffman_str

def combineHuffman(huffman_str_list, idnex_len, add_len):
    index_huffman_dic = {i[:idnex_len]:i[idnex_len:] for i in huffman_str_list}
    idnex_list = list(index_huffman_dic.keys())
    idnex_list.sort()

    info_huffman_list = [index_huffman_dic[i] for i in idnex_list]


    total_huffman_str = info_huffman_list[0]
    len_info = int(len(total_huffman_str)/4)
    for info_huffman in tqdm(info_huffman_list[1:]):
        total_huffman_str += info_huffman[-len_info:]

    if add_len != 0:
        total_huffman_str = total_huffman_str[:-add_len]

    return total_huffman_str 

def huffmanToByte(total_huffman_str):
    byte_list = []
    start = 0
    n=0
    while True:
        hufffman_code = total_huffman_str[start: start+5]
        step = 5
        if hufffman_code > "22201":
            hufffman_code = total_huffman_str[start: start+6]
            step = 6

        if len(hufffman_code) <5:
            break

        start += step
        decimal_num =int(hufffman_code, 3)
        if step == 6:
            decimal_num -= 472
        byte_list.append(decimal_num)
        # print(n, step)
        # if n == 51:
        #     print(hufffman_code, step, decimal_num)
        #     break
        n += 1
    return byte_list
        # _bin_str = bin(decimal_num)[2:]
        # _bin_str = "0"*(8-len(_bin_str)) + _bin_str
        # bin_str

def saveResult(byte_list, save_path):
    # print(byte_list)
    with open(save_path, "wb") as f:
        for i in byte_list:
            f.write(bytes([i]))


def goldmanDecode(nt_list, save_path, idnex_len, add_len):
    huffman_str_list = []
    for nt_seq in tqdm(nt_list):
        huffman_str = decodeNt(nt_seq)
        huffman_str_list.append(huffman_str)

    total_huffman_str = combineHuffman(huffman_str_list, idnex_len, add_len)
    # return total_huffman_str
    byte_list = huffmanToByte(total_huffman_str)
    # return byte_list
    # print("byte_list", len(byte_list)) #bkp
    saveResult(byte_list, save_path)
    

def readInput(input_path):
    nt_list = []
    with open(input_path) as f:
        for i in f:
            nt_seq = i.strip()
            nt_list.append(nt_seq)
    return nt_list

def goldmanDecodeMain(input_path, save_path, idnex_len, add_len):
    nt_list = readInput(input_path)
    decode_r = goldmanDecode(nt_list, save_path, idnex_len, add_len)
    return decode_r
