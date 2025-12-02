from scipy.constants import Avogadro
from itertools import groupby
import random
import importlib
import numpy as np 

from dnabyte.oligo import Oligo, complement, translate_nested_list, translate_nested_list_poly,back_translate_nested_list_poly_binom_real, back_translate_nested_list_chain, back_translate_nested_list_poly, back_translate_nested_list, back_translate_nested_list_poly_binom,back_translate_nested_list_real,back_translate_nested_list_chain_real
from dnabyte.oligopool import OligoPool
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
import numpy as np

class SimulateSynthesis:
    """
    Simulate the synthesis of single stranded oligos into long chaines of double stranded DNA.
    """
    def __init__(self, params, logger=None):
            
            self.assembly_structure = params.assembly_structure
            self.encoding_method = params.encoding_method
            self.synthesis_method = params.synthesis_method
            self.params = params
            self.synthesis_plugins = params.synthesis_plugins

            # To allow simulation with large quantites of oligos, we will scale down the number, perform
            # the simulation, and then scale back up the results.
            # if self.synthesis_method != None:
            #     scale = 1
            #     while self.mean > 1000:
            #         self.mean = self.mean / 10
            #         self.std_dev = self.std_dev / 10
            #         scale = scale * 10

            #     self.scale = scale

    def simulate(self, data):
        """
        Simulate the assembly of DNA sequences.
        
        :param data: A list of DNA sequences.
        :return: A list of assembled DNA sequences.
        """
        
        if self.assembly_structure == 'synthesis':
            if self.synthesis_method == None:
                assembled_data = InSilicoDNA(data.data)
                info = {"copy_number": 0}
                return assembled_data, info
            else:
                if isinstance(data, NucleobaseCode):
                    # Dynamically find the appropriate class based on synthesis_method
                    try:
                        if self.params.synthesis_method in [3, 4, 5, 6, 7, 68, 69, 70, 71]:
                            synthesis_class = self.synthesis_plugins["mesa"]
                            plugin = synthesis_class(self)
                        else:
                            synthesis_class = self.synthesis_plugins[self.synthesis_method.lower()]
                            plugin = synthesis_class(self)  # Instantiate the plugin class

                        # Call the simulate method of the plugin class
                        data_seq, info = plugin.simulate(data.data)
                        data_seq = InSilicoDNA(data_seq)
                        # TODO: change the initialization of data_seq to be more elegant
                        #data_seq.data = data_seq
                        return data_seq, info
                        
                    except KeyError:
                        raise ValueError(f"Synthesis method '{self.synthesis_method}' not recognized.")
                else:
                    raise ValueError("The input data is not an instance of EncodedData.")
        else:
            print(data)
            if isinstance(data, NucleobaseCode):
                assembly = importlib.import_module(f"dnabyte.encoding.{self.encoding_method}.assembly")
                assembled_data, info = assembly.assembly(data, self.params)
                assembled_data = InSilicoDNA(assembled_data)
                return assembled_data, info
        

    def print_list_structure(self, lst, level=0):
        """
        Recursively prints the structure of a list of lists.

        Args:
            lst (list): The list to analyze.
            level (int): The current nesting level (used for indentation).
        """
        if isinstance(lst, list):
            print("  " * level + f"Level {level}: List with {len(lst)} elements")
            for item in lst:
                self.print_list_structure(item, level + 1)
        else:
            if isinstance(lst, Oligo):
                print("  " * level + f"Level {level}: {type(lst).__name__} ({lst})")
            elif isinstance(lst, OligoPool):
                for oligo in lst.pool:
                    print("  " * level + f"Level {level}: {type(oligo).__name__} ({oligo})")
                
    def process_tuple_list(tuple_list):
        def flip_tuple(t):
            return t[::-1]

        # Ensure the tuple containing 'A' is the first tuple
        if 'A' == tuple_list[1][0] or 'A' == tuple_list[1][-1]:
            tuple_list = (tuple_list[1], tuple_list[0])

        # Ensure 'A' is the first entry of the first tuple
        if tuple_list[0][0] != 'A':
            tuple_list = (flip_tuple(tuple_list[0]), flip_tuple(tuple_list[1]))

        return tuple_list
    
    def add_empty_to_innermost(self, lst, library):
        finishedlist = []
        finishedlist.append(Oligo([complement(lst[0].motifs[0],library.dictmotives), 'empty']))
        finishedlist.extend(lst)
        return finishedlist

    def print_nested_list(self,nested_list, level=0):
        for item in nested_list:
            if isinstance(item, list):
                self.print_nested_list(item, level + 1)
            else:
                print(item,  'oligo')

    def nest_list(self, lst):
        """
        Recursively nests a list into the structure [[[[a, b], c], d], e].

        Args:
            lst (list): The input list to be nested.

        Returns:
            list: The nested list structure.
        """
        if len(lst) <= 1:
            return lst  # Base case: if the list has 1 or fewer elements, return it as is
        return [self.nest_list(lst[:-1]), lst[-1]]  # Recursively nest all but the last element with the last element

    def _validate(self, params):
        """
        Validates the provided parameters for encoding and assembly.

        This method checks if the provided parameters for encoding and assembly
        are valid according to predefined criteria.

        Parameters:
        params (Params): The parameters to be validated.

        Raises:
        ValueError: If any of the parameters are invalid.
        """
        # Implement validation logic here
        pass

