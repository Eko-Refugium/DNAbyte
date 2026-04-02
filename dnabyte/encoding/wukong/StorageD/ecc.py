# -*- coding: utf-8 -*-
from reedsolo import RSCodec, ReedSolomonError
from tqdm import tqdm
import math
import logging

log = logging.getLogger('mylog')

bin_to_str_list = list()
for i in range(256):
    bin_to_str_list.append(bin(i)[2:].rjust(8,'0'))

class ReedSolomon():
    def __init__(self, check_bytes:int):
        """
        reed-solomon code
        @param check_bytes: the number of rs code (bytes)
        """
        self.check_bytes = check_bytes
        self.tool = RSCodec(check_bytes)
        self.remainder = 0
        log.debug('ecc check_bytes:{}'.format(check_bytes))
        
    def insert_one(self, binstring, group=1):
        # If the binary length is not a multiple of 8, it will be left-padded with zeros
        self.remainder = 0
        if len(binstring)%8!=0:
            self.remainder = 8 - len(binstring) % 8
            binstring = '0'*self.remainder + binstring

        output_string = ''
        if group>=1 and type(group)==int:
            if group<=len(binstring)/8/255:
                group = 1
            group_len = math.ceil(len(binstring)/8/group)*8
            for i in range(group):
                last = group_len*(i+1)
                if last>=len(binstring): last = None
                group_string = binstring[group_len*i:last]
                byte_list = bytearray([int(group_string[i:i + 8], 2) for i in range(0, len(group_string), 8)])
                rs_list = self.tool.encode(byte_list)
                str_list = [bin_to_str_list[i] for i in rs_list]
                output_string += ''.join(str_list)
        else:
            raise Exception("Wrong rs group")
        return output_string[self.remainder:]
    
    
    def remove_one(self, binstring, group=1):
        self.remainder = 0
        # If the binary length is not a multiple of 8, it will be left-padded with zeros
        if len(binstring)%8!=0:
            self.remainder = 8 - len(binstring) % 8
            binstring = '0'*self.remainder + binstring

        data_len = len(binstring) - group*self.check_bytes*8
        output_string = ""
        if group>=1 and type(group)==int:
            if group<=data_len/8/255:
                group = 1
            group_len = math.ceil(data_len/8/group)*8 + self.check_bytes*8
            for i in range(group):
                last = group_len*(i+1)
                if last >= len(binstring): last = None
                group_string = binstring[i*group_len:last]
                byte_list = bytearray([int(group_string[i:i + 8], 2) for i in range(0, len(group_string), 8)])
                try:
                    decode_byte_list, full_list, err_pos = self.tool.decode(byte_list)
                    if len(err_pos)!=0:
                        print(err_pos)
                    str_list = [bin_to_str_list[i] for i in decode_byte_list]
                    output_string += ''.join(str_list)
                except ReedSolomonError:
                    # Irreparable
                    return {"data": binstring, "type": False}
                except IndexError:
                    # No data acquisition
                    return {"data": binstring, "type": False}
            return {"data": output_string, "type": True}
        else:
            raise Exception("Wrong rs group")
    
    
    def insert(self, segment_list, group=1):
        log.debug('start adding rs code.')
        res_list = list()
        if len(segment_list)==0:
            raise ValueError("Empty data.")
        # Insert rs codes into multiple sequences
        if type(segment_list)==list and type(segment_list[0]==str):
            pro_bar = tqdm(total=len(segment_list), desc="Add RSCode")
            for segment in segment_list:
                res = self.insert_one(segment, group=group)
                res_list.append(res)
                pro_bar.update()
            pro_bar.close()
        else:
            raise Exception("Input error")
        return res_list
    
    
    def remove(self, segment_list, oriLen, group=1):
        log.debug('Check and  remove rs code.')
        if len(segment_list)==0:
            raise ValueError("Empty data.")
        bit_segments = []
        error_bit_segments = []
        error_indices = []
        if type(segment_list) == list and type(segment_list[0]) == str:
            error_rate = 0
            pro_bar = tqdm(total=len(segment_list), desc="Del RSCode")
            for index, segment in enumerate(segment_list):
                if segment is not None:
                    output = self.remove_one(segment, group=group)
                    data, is_succeed = output.get("data"), output.get("type")
                    if is_succeed and len(data)>=oriLen:
                        bit_segments.append(data[len(data) - oriLen:])
                    else:
                        error_rate += 1
                        error_indices.append(index)
                        error_bit_segments.append(data)
                else:
                    error_rate += 1
                    error_indices.append(index)
                    error_bit_segments.append(None)
                pro_bar.update()
            pro_bar.close()
            error_rate /= len(segment_list)
        else:
            raise ValueError("Input error")
        return {"bit": bit_segments, "e_r": error_rate, "e_i": error_indices, "e_bit": error_bit_segments}
