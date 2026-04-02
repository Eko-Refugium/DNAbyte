# -*- coding: utf-8 -*-
from abc import ABC,abstractmethod
from os import path,stat
from typing import Tuple,List
from datetime import datetime
import logging
log = logging.getLogger('mylog')

from dnabyte.encoding.wukong.StorageD.tools import  EncodeParameter, DecodeParameter, FileTools, SplitTools, RsTools, BaseTools, CodecException
from dnabyte.encoding.wukong.StorageD.getPrimerPair import PrimerDesign, getOther

base_list = ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']

class AbstractEncode(ABC):
    """
    Abstract encode class
    """
    def __init__(self, input_file_path:str, output_dir:str, tool_name,seq_bit_to_base_ratio=1, index_redundancy=0, **encode_parameters):
        """
        Subclass may add some new parameters, 
        and add 'super().__init__(input_file_path, output_dir, **kwargs)' in its __init__.
        """
        self.tool_name = tool_name # in filename
        self.input_file_path = input_file_path
        self.output_dir = output_dir if output_dir[-1] in ['/','\\'] else output_dir+'/'
        if not path.exists(self.output_dir):
            raise CodecException('Output dir not exist : {}'.format(self.output_dir))
        self.file_base_name,self.file_extension =  FileTools.get_file_name_extension(self.input_file_path)
        self.output_file_path = self.output_dir + self.file_base_name + "_{}.fasta".format(tool_name)
        self.codec_param=EncodeParameter(**encode_parameters)
        self.seq_bit_to_base_ratio = seq_bit_to_base_ratio
        self.index_redundancy = index_redundancy
        self.total_bit =  8*stat(self.input_file_path).st_size
        self.seq_num = 0
        self.bin_split_len = 0
        self.index_length = 0
        self.run_time = '0:00:00'
        self.encode_time = '0:00:00'
        self.density = 0
        self.total_base = 0
        self.rs_group = 0
        self._check_params()
        self.codec_worker = None

    def _check_params(self):
        if self.codec_param.max_homopolymer<=0:
            raise CodecException("homopolymer:{} must be greater than 0".format(self.codec_param.homopolymer))
        if self.codec_param.min_gc<=0 or self.codec_param.max_gc>=1:
            raise CodecException("min_gc:{} must be greater than 0 and max_gc:{} must be less than 1".format(self.codec_param.min_gc, self.codec_param.max_gc))
        if self.codec_param.min_gc > self.codec_param.max_gc:
            raise CodecException("min_gc:{} could not be greater than max_gc:{}".format(self.codec_param.min_gc, self.codec_param.max_gc))
        if self.codec_param.primer_length<18 and self.codec_param.primer_length>24: # set in primer3Settiong.py
            raise CodecException("primer_length:{} should be between 18 and 24")
        if self.codec_param.sequence_length<=0:
            raise CodecException("sequence_length must be greater than 0")
        
    
    def _set_density(self):
        if self.total_base!=0:
            self.density = self.total_bit/self.total_base
        else:
            log.error("set density error, total base is 0")
    
    @abstractmethod
    def _encode(self, bit_segments: list) -> list:
        """
        This part is implemented according to different algorithms
        Processed Binary Sequence(str) List >> Base Sequence(str) List
        """
        # return base_segments
        pass
    
    def _file_to_bin(self):
        """
        get binary string from file
        """
        bin_str = FileTools.file_to_bin(self.input_file_path)
        return bin_str

    def _split_bin(self, binstring:str) -> Tuple[int,int,List]:
        """
        split binstring to bin segments
        """
        index_length, bin_split_len, seq_num, rs_group = SplitTools.get_bin_split_length(len(binstring),
                                                                               self.seq_bit_to_base_ratio,
                                                                               self.codec_param,
                                                                               index_redundancy=self.index_redundancy)
        bit_segments = SplitTools.split(binstring, bin_split_len, self.codec_param.add_redundancy,
                                        index_redundancy=self.index_redundancy)
        return bin_split_len, index_length, rs_group, bit_segments
        
    def _add_rscode(self, bit_segments):
        """
        add rscode
        """
        log.info('add rs')
        res_segments = RsTools.add_rs(bit_segments, self.codec_param.rs_num, self.rs_group)
        return res_segments
        
    def _add_primer(self):
        log.info('add primer')
        primer_designer = PrimerDesign([self.output_file_path], iBreakNum=1, iBreakSec=300,iPrimerLen=self.codec_param.primer_length, sBaseDir=self.output_dir, bTimeTempName=False, bLenStrict=True)
        res = primer_designer.getPrimer()
        if len(res)>0:
            left_primer = res.left[0]
            right_primer = res.right[0]
            right_primer_reverse = getOther(right_primer) # add reverse complementary at the right
            pool = []
            with open(self.output_file_path, 'r') as file:
                param = self._set_param_line(left_primer, right_primer)
                pool.append(param)
                index = 0
                for line in file.readlines()[1:]:
                    if line[0] in base_list:
                        pool.append(">seq_{}\n".format(index+1) + left_primer + line.strip('\n') + right_primer_reverse + '\n')
                        index += 1
            with open(self.output_file_path, 'w') as f:
                f.writelines(pool)
            return self.output_file_path
        else:
            log.error("No primer")
            return self.output_file_path
    
    def _set_param_line(self, left_primer="", right_primer=""):
        param = ">totalBit:{},binSegLen:{},leftPrimer:{},rightPrimer:{},fileExtension:{},bRedundancy:{},RSNum:{}\n".format(
            self.total_bit, self.bin_split_len, left_primer, right_primer, self.file_extension,
            int(self.codec_param.add_redundancy), int(self.codec_param.rs_num))
        return param
        
    def common_encode(self):
        """
        common encode flow
        """
        self._check_params()
        tm_run = datetime.now()
        log.debug('read file')
        bin_str = self._file_to_bin()
        # split
        self.bin_split_len, self.index_length, self.rs_group, bit_segments = self._split_bin(bin_str)
        self.seq_num = len(bit_segments)
        #add rscode
        if self.codec_param.rs_num>0:
            bit_segments = self._add_rscode(bit_segments)
            
        # encode
        tm_encode = datetime.now()
        base_segments = self._encode(bit_segments)
        self.encode_time = str(datetime.now() - tm_encode)
        self.total_base = len(base_segments)*len(base_segments[0])
        
        param = self._set_param_line()
        log.debug('write')
        with open(self.output_file_path, 'w') as f:
            f.write(param)
            for index,seg in enumerate(base_segments):
                f.write(">seq_{}\n".format(index+1) + seg + '\n')
        
        # add primer (blast need fasta file to compare)
        if self.codec_param.add_primer:
            self._add_primer()
            self.total_base += (len(base_segments)*2*self.codec_param.primer_length)
            
        self.run_time = str(datetime.now()-tm_run)
        self._set_density()
        return self.output_file_path
    
    
class AbstractDecode(ABC):
    def __init__(self, input_file_path:str, output_dir:str, index_redundancy=0):
        """
        Subclass may add some new parameters, 
        and add 'super().__init__(input_file_path, output_dir)' in its __init__.
        """
        self.input_file_path = input_file_path
        self.output_dir = output_dir if output_dir[-1] in ['/','\\'] else output_dir+'/'
        if not path.exists(self.output_dir):
            raise CodecException('Output dir not exist')
        self.file_base_name,self.file_extension =  path.splitext(path.basename(self.input_file_path))
        self.rs_err_rate = 0.0 # error correction failed
        self.rs_err_indexs = [] # error correction failed
        self.repaired_rate = 0.0 # missing but repaired segments
        self.miss_err_rate = 0.0 # missing and unrepaired segments
        self.repaired_indexs = [] # missing but repaired segments
        self.miss_err_indexs = [] # missing and unrepaired segments
        self.run_time = 0.0
        self.index_redundancy = index_redundancy
        self.codec_worker = None
    
    def _check_file_param(self, param_list:list):
        with open(self.input_file_path,'r') as file:
            first_line = file.readline().strip()
            param_line = first_line.strip(">").strip(';').split(',')

        param_dict = dict()
        for param in param_line:
            temp = param.split(":")
            param_dict[temp[0]] = temp[1]
        list_input = set(param_dict.keys())
        list_default = set(param_list)
        err_list = list(list_input.difference(list_default)) + list(list_default.difference(list_input))
        log.info(err_list)
        if len(err_list)!=0:
            log.error("wrong para : {}".format((',').join(err_list)))
            raise CodecException("Failed to parse parameter. Please Check that the input file and the selected codec algorithm are consistent.")
        return param_dict

    def _parse_param(self):
        """
        parse parameter from init input and input_file
        """
        param_dict = self._check_file_param(["totalBit","binSegLen","leftPrimer","rightPrimer","fileExtension","bRedundancy","RSNum"])  
        param = {
            'total_bit': int(param_dict['totalBit']),
            'bin_seg_len' : int(param_dict['binSegLen']),
            'add_redundancy' : bool(int(param_dict['bRedundancy'])),
            'rs_num' :int(param_dict['RSNum']),
            'left_primer_len' : len(param_dict['leftPrimer']),
            'right_primer_len' : len(param_dict['rightPrimer']),
            'file_extension' : param_dict['fileExtension']
            }
        self.codec_param = DecodeParameter(**param)
        self.index_length,self.seq_num = self._get_indexlen_splitnum()
        self.ori_len = self.codec_param.bin_seg_len + self.index_length
        self.output_file_path = self.output_dir + self.file_base_name + "_decode" + self.codec_param.file_extension

    def _get_indexlen_splitnum(self):
        index_length, seq_num = SplitTools.get_indexlen_and_segnum(self.codec_param.total_bit,
                                                                     self.codec_param.bin_seg_len,
                                                                     self.codec_param.add_redundancy,
                                                                     index_redundancy=self.index_redundancy)
        return index_length, seq_num
    
    def _get_base_line_list(self):
        base_line_list = []
        with open(self.input_file_path,'r') as file:
            if self.codec_param.right_primer_len==0:
                def get_line():
                    base_line_list.append(line.strip()[self.codec_param.left_primer_len:])
            else:
                def get_line():
                    base_line_list.append(line.strip()[self.codec_param.left_primer_len:-self.codec_param.right_primer_len])
            # del primer
            for line in file.readlines():
                if line[0] in base_list:
                    get_line()
        base_line_list = [x.strip() for x in base_line_list]
        return base_line_list
        
    @abstractmethod
    def _decode(self, base_line_list:list()) -> list():
        """
        The input parameter is a list of base sequences, and the return value is a list of binary sequences
        """
        # return bit_segments
        pass
    
    def _del_rscode(self, bit_segs):
        log.debug('del rscode')
        res_dict = RsTools.del_rs(bit_segs, self.ori_len, self.codec_param.rs_num)
        err_segs = [bit_segs[i] for i in res_dict.get("err_index")]
        self.rs_err_rate = res_dict.get("err_rate")
        err_lists_tmp = [int(bit_segs[i][:self.index_length],2) for i in res_dict["err_index"]]
        for index in err_lists_tmp:
            if index < self.seq_num:
                self.rs_err_indexs.append(index)
        return res_dict['validate_bit_seg'], err_segs

    def _repair_segment(self, bit_segs):
        res_dict = SplitTools.repair_segment(self.index_length, bit_segs, self.seq_num,
                                             self.codec_param.add_redundancy)
        self.repaired_rate = res_dict['repaired_rate']
        self.repaired_indexs = res_dict['repaired_indexs']
        self.miss_err_rate = res_dict['failed_rate']
        self.miss_err_indexs = res_dict['failed_indexs']
        return res_dict["res_segments"]
    
    def _bin_to_file(self, bit_str):
        FileTools.bin_to_file(bit_str, self.output_file_path)
    
    def common_decode(self):
        self._parse_param()

        tm_run = datetime.now()
        base_line_list = self._get_base_line_list()
        tm_decode = datetime.now()
        bit_segments = self._decode(base_line_list)
        self.decode_time = str(datetime.now() - tm_decode)
        log.debug('sort')
        sorted_binstr = BaseTools.sort_segment(bit_segments, self.index_length)
        sorted_binstr = sorted_binstr[:self.seq_num]

        # del rscode
        if self.codec_param.rs_num>0:
            validate_bit_segs, err_bit_segs = self._del_rscode(sorted_binstr)
        else:
            validate_bit_segs, err_bit_segs = sorted_binstr, []
        # TODO: err_bit_segs and self.rs_err_indexs may be used to repaire
        #repair and count missing segments
        log.debug('repair')
        validate_bit_segs = self._repair_segment(validate_bit_segs)
        
        log.debug('merge')
        res_bit_str = SplitTools.merge(validate_bit_segs)
        res_bit_str = res_bit_str[:self.codec_param.total_bit]
        
        log.debug('write')
        self._bin_to_file(res_bit_str)
        
        self.run_time = str(datetime.now() - tm_run)
        return self.output_file_path