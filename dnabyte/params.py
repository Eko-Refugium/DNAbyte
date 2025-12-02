import os
import importlib
import inspect

from dnabyte.encode import Encode
from dnabyte.store import SimulateStorage
from dnabyte.load_plugins import load_plugins

from dnabyte.binarize import Binarize
from dnabyte.sequence import SimulateSequencing
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.library import Library
from typing import List, Tuple
from pyxdameraulevenshtein import damerau_levenshtein_distance

class Params:

    def __init__(self, debug=False, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.debug = debug
        print(self)
        # Load plugins
        self.encoding_plugins, self.synthesis_plugins, self.storage_plugins, self.sequencing_plugins, self.binarization_plugins = load_plugins(self.encoding_method, self.synthesis_method, self.storage_conditions, self.sequencing_method, self.binarization_method)


        print(self.encoding_method,self.encoding_plugins, 'encoding method in params')
        
        if self.encoding_method is None:
            raise ValueError("Encoding method must be specified")
        elif self.encoding_method in self.encoding_plugins:
            encoding = importlib.import_module(f"dnabyte.encoding.{self.encoding_method}.encode")
            attributes_enc = encoding.attributes(self) 
            for keys, value in attributes_enc.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid encoding method: {self.encoding_method}")
        self.encoding_plugins, self.synthesis_plugins, self.storage_plugins, self.sequencing_plugins, self.binarization_plugins = load_plugins(self.encoding_method, self.synthesis_method, self.storage_conditions, self.sequencing_method, self.binarization_method)
        if self.synthesis_method is None:
            pass
        elif 'synthesis_' + self.synthesis_method in self.synthesis_plugins:
            synthesis = importlib.import_module(f"dnabyte.synthesis.synthesis_{self.synthesis_method}")
            attributes_synth = synthesis.attributes(self)
            for keys, value in attributes_synth.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid synthesis method: {self.synthesis_method}")
        print(self.storage_conditions,self.storage_plugins, 'storage method in params')
        if self.storage_conditions is None:
            pass
        elif 'storage_' + self.storage_conditions in self.storage_plugins:
            storage = importlib.import_module(f"dnabyte.storage.storage_{self.storage_conditions}")
            attributes_store = storage.attributes(self)
            for keys, value in attributes_store.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid storage conditions: {self.storage_conditions}")
        
        if self.sequencing_method is None:
            pass
        elif 'sequencing_' + self.sequencing_method in self.sequencing_plugins:
            sequencing = importlib.import_module(f"dnabyte.sequencing.sequencing_{self.sequencing_method}")
            attributes_seq = sequencing.attributes(self)
            for keys, value in attributes_seq.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid sequencing method: {self.sequencing_method}")
        
        print(self.binarization_method,self.binarization_plugins, 'binarization method in params')
        if self.binarization_method is None:
            pass
        elif self.binarization_method in self.binarization_plugins:
            binarization = importlib.import_module(f"dnabyte.binarization.{self.binarization_method}")
            attributes_bin = binarization.attributes(self)
            for keys, value in attributes_bin.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid binarization method: {self.binarization_method}")
        
        print(self)
        # breakpoint()
        

        # if self.inner_error_correction not in ['ltcode', None]:
        #     raise ValueError("Invalid inner error correction")
        
        # if self.outer_error_correction not in ['reedsolomon', None]:
        #     raise ValueError("Invalid outer error correction")

        # if self.outer_error_correction == 'reedsolomon':
        #     if self.reed_solo_percentage < 0 or self.reed_solo_percentage > 100:
        #         raise ValueError("Reedsolomon percentage not compatible")
            
        # if self.inner_error_correction == 'ltcode':
        #     if self.percent_of_symbols < 0 or self.percent_of_symbols > 100:
        #         raise ValueError("Percent of symbols not compatible")
            
        #     if self.index_carry_length < 0:
        #         raise ValueError("Index carry length not compatible")    
        
        # if self.assembly_structure == 'positional_assembly' or self.assembly_structure == 'linear_assembly':
        #     if self.library is None or self.library_name is None:
        #         raise ValueError("Library not provided")
            
        # if self.encoding_scheme == 'binomial_encoding':
        #     if self.sigma_amount is None or self.theory is None:
        #         raise ValueError("Missing parameters for encoding scheme")
            
        # if self.encoding_scheme == 'linear_encoding' and self.assembly_structure == 'linear_assembly':
        #     if self.codeword_length is None or self.theory is None:
        #         raise ValueError("Missing parameters for encoding scheme")
            
        # if (self.assembly_structure == 'synthesis' or (self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'linear_encoding')):
        #     if self.codeword_length is None:
        #         raise ValueError("Missing parameters for encoding scheme")
            
        # if self.encoding_scheme == 'linear_encoding' and self.assembly_structure == 'positional_assembly':
        #     if self.theory is None:
        #         raise ValueError("Missing parameters for encoding scheme")

        # if self.sequencing_method not in [41, 40, 37, 36, 39, 38, 35,'iid', None]:
        #     raise ValueError("Invalid sequencing method")
        
        # if self.synthesis_method not in [3, 4, 5, 6, 7, 68, 69, 70, 71, 'nosynthpoly', None]:
        #     raise ValueError("Invalid Synthesis Error")
        
        # # Validate storage conditions
        # if isinstance(self.storage_conditions, list):
        #     for condition in self.storage_conditions:
        #         if condition not in self.storage_plugins and condition not in ['permafrost', 'roomtemperature', 'biogene', None]:
        #             raise ValueError("Invalid storage conditions")
        # elif self.storage_conditions is None:
        #     pass
        # elif isinstance(self.storage_conditions, str):
        #     if self.storage_conditions not in self.storage_plugins and self.storage_conditions not in ['permafrost', 'roomtemperature', 'biogene', None]:
        #         raise ValueError("Invalid storage conditions")
        # else:
        #     raise ValueError("Invalid storage conditions")

        
    def __str__(self):
        output = "Params:\n"
        for key, value in self.__dict__.items():
            output += f"  {key}: {value}\n"
        return output


    @classmethod
    def params_range(cls, **kwargs):
        """
        Alternative constructor that creates a list of Params instances for each item in the list if one of the 
        parameters is a non-empty list.
        
        Parameters:
        kwargs (dict): The parameters for the Params instances.
        
        Returns:
        list: A list of Params instances.
        """
        params_list = []
        
        # Find the first attribute that is a non-empty list
        for attr_name, attr_value in kwargs.items():
            if isinstance(attr_value, list) and attr_value:
                for item in attr_value:
                    new_kwargs = {**kwargs, attr_name: item}
                    params_list.append(cls(**new_kwargs))
                break
        else:
            # If no attribute is a non-empty list, return a list with a single instance
            params_list.append(cls(**kwargs))
        
        return params_list

