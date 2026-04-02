# -*- coding: utf-8 -*-
from datetime import datetime
from math import ceil
import json

from dnabyte.encoding.wukong.StorageD.abstract_codec import AbstractEncode,AbstractDecode
from dnabyte.encoding.wukong.StorageD.tools import  DecodeParameter, CodecException
from dnabyte.encoding.wukong.StorageD.wukong import Wukong
from dnabyte.encoding.wukong.StorageD.church import churchEncode
from dnabyte.encoding.wukong.StorageD.churchDecode import churchDecode
from dnabyte.encoding.wukong.StorageD.goldman import goldmanEncode
from dnabyte.encoding.wukong.StorageD.goldmanDecode import goldmanDecode

import logging
log = logging.getLogger('mylog')

class WukongEncode(AbstractEncode):
    def __init__(self, input_file_path:str, output_dir:str, sequence_length:int, rule_num=0,
                 max_homopolymer=6, min_gc=0.3, max_gc=0.7, rs_num=0,  add_redundancy=False, 
                 add_primer=False, primer_length=20, same_res=True):
        self.rule_num = rule_num
        self.same_res = same_res
        self.virtual_segment = 0
        index_redundancy = 1 #for virtual segments in wukong
        seq_bit_to_base_ratio = 1
        data_len = sequence_length-primer_length if add_primer else sequence_length
        if data_len%2!=0:
            sequence_length -= 1
            raise CodecException("Wukong: The sequence length after removing the Flanking Sequence should be an even number")
        super().__init__(input_file_path, output_dir, 'wukong', seq_bit_to_base_ratio=seq_bit_to_base_ratio,
                         index_redundancy=index_redundancy,sequence_length=sequence_length, max_homopolymer=max_homopolymer, 
                         min_gc=min_gc, max_gc=max_gc, rs_num=rs_num, add_redundancy=add_redundancy, add_primer=add_primer, 
                         primer_length=primer_length)

    def _encode(self, bit_segments):
        self.codec_worker = Wukong(rule_num=self.rule_num)
        base_segments = self.codec_worker.wukong_encode(index_length=self.index_length, bit_segments=bit_segments, max_homopolymer=self.codec_param.max_homopolymer,
                                                    max_content=self.codec_param.max_gc, min_content=self.codec_param.min_gc, hasseed=self.same_res)
        self.virtual_segment = self.codec_worker.addtion_num
        return base_segments
    
class WukongDecode(AbstractDecode):
    def __init__(self, input_file_path:str, output_dir:str, rule_num=0):
        self.rule_num = rule_num
        index_redundancy = 1
        super().__init__(input_file_path, output_dir, index_redundancy=index_redundancy)

    def _decode(self, base_line_list):
        self.codec_worker = Wukong(rule_num=self.rule_num)
        bit_segments = self.codec_worker.wukong_decode(base_line_list)
        return bit_segments

############################################  
"""
Reference:
* Church, G. M., et al. (2012). "Next-generation digital information storage in DNA." Science 337(6102): 1628.
* DOI: 10.1126/science.1226355
* https://pubmed.ncbi.nlm.nih.gov/22903519/
""" 
class ChurchEncode(AbstractEncode):
    def __init__(self, input_file_path:str, output_dir:str, sequence_length:int, max_homopolymer=6, 
                 rs_num=0, add_redundancy=False, add_primer=False, primer_length=20):        
        index_redundancy = 0 #for virtual segments in wukong
        seq_bit_to_base_ratio = 1
        super().__init__(input_file_path, output_dir, 'church', seq_bit_to_base_ratio=seq_bit_to_base_ratio, index_redundancy=index_redundancy,
                         sequence_length=sequence_length, max_homopolymer=max_homopolymer, rs_num=rs_num, add_redundancy=add_redundancy, 
                         add_primer=add_primer, primer_length=primer_length)
    
    def _encode(self, bit_segments):
        base_segments = churchEncode(bit_segments, rep_num=self.codec_param.max_homopolymer)
        return base_segments
    

class ChurchDecode(AbstractDecode):
    def __init__(self, input_file_path:str, output_dir:str):
        super().__init__(input_file_path, output_dir)

    def _decode(self, base_line_list):
        bit_segments = churchDecode(base_line_list)
        return bit_segments
    
############################################ 
"""
Reference:
* Goldman, N., et al. (2013). "Towards practical, high-capacity, low-maintenance information storage in synthesized DNA." Nature 494(7435): 77-80.
* DOI: 10.1038/nature11875
* https://www.ncbi.nlm.nih.gov/pubmed/23354052
"""
class GoldmanEncode(AbstractEncode): 
    def __init__(self, input_file_path:str, output_dir:str, sequence_length:int):
        """
        Subclass may add some new parameters, 
        and add 'super().__init__(input_file_path, output_dir, **kwargs)' in its __init__.
        """
        super().__init__(input_file_path, output_dir, 'goldman', **{'sequence_length':sequence_length})
        
    
    def _encode(self, bit_segments: str) -> list:
        pass

    # rewrite
    def common_encode(self):
        tm_run = datetime.now()
        log.debug('read')
        bin_str = self._file_to_bin()
        # encode
        log.debug('encode')
        tm_encode = datetime.now()
        nt_seq_list, self.index_length, self.add_len, ternary_seg_list, ternary_str = goldmanEncode(bin_str, self.codec_param.sequence_length)
        self.encode_time = str(datetime.now()-tm_encode)
        
        self.total_base = len(nt_seq_list)*len(nt_seq_list[0])
        self.seq_num = len(nt_seq_list)
        param = ">indexLen:{},addLen:{},fileExtension:{}\n".format(
            self.index_length, self.add_len, self.file_extension)
        with open(self.output_file_path, 'w') as f:
            f.write(param)
            for index,seg in enumerate(nt_seq_list):
                f.write(">seq_{}\n".format(index+1) + seg + '\n')
        
        self.run_time = str(datetime.now()-tm_run)
        self._set_density()
        return self.output_file_path
        

class GoldmanDecode(AbstractDecode):
    def __init__(self, input_file_path:str, output_dir:str):
         super().__init__(input_file_path, output_dir)
         
    def _decode(self):
        pass

    # rewrite
    def _parse_param(self):
        param_dict = self._check_file_param(["indexLen","addLen","fileExtension"])
        self.index_length = int(param_dict['indexLen'])
        self.add_len = int(param_dict['addLen'])
        self.file_extension = param_dict['fileExtension']
        self.codec_param=DecodeParameter(file_extension=self.file_extension)
        self.output_file_path = self.output_dir + self.file_base_name + "_decode" + self.codec_param.file_extension
        
    def common_decode(self):
        self._parse_param()
        tm_run = datetime.now()
        base_line_list = self._get_base_line_list()
        #decode
        goldmanDecode(base_line_list, self.output_file_path, self.index_length, self.add_len)
        self.run_time = str(datetime.now() - tm_run)
        return self.output_file_path