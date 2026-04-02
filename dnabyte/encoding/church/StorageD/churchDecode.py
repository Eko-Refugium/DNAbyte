# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 10:37:21 2021

@author: Wei QWiang
"""

import argparse
from tqdm import tqdm

map_dic = {"A":"0", "C":"0", "G":"1", "T":"1"}

def linuxCommand():
    parser = argparse.ArgumentParser(description = "Church's algorithm in the pipeline of ATOM.")
    parser.add_argument('i')
    parser.add_argument('o')
    args = parser.parse_args()

    return args.i, args.o


def seqDecode(nt_seq: str, map_dic : dict = map_dic) -> str:
    bin_str = ""
    for nt in nt_seq:
        bin_str += map_dic[nt]

    return bin_str

def churchDecode(seq_list: list, map_dic : dict = map_dic) -> list:
    bin_list = []
    pro_bar = tqdm(total=len(seq_list), desc="Decoding")
    for nt_seq in seq_list:
        bin_str = seqDecode(nt_seq)
        bin_list.append(bin_str)
        pro_bar.update()
    pro_bar.close()

    return bin_list

def churchDecodeMain(input_path, output_path, map_dic : dict = map_dic):
    o_f = open(output_path, "w")
    with open(input_path) as f:
        for i in f:
            nt_seq = i.strip()
            bin_str = seqDecode(nt_seq)
            o_f.write(bin_str+"\n")

    o_f.close()

