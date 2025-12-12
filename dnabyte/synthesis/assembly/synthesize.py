
import importlib

from dnabyte.synthesize import SimulateSynthesis
from .oligo import Oligo, complement, translate_nested_list, translate_nested_list_poly,back_translate_nested_list_poly_binom_real, back_translate_nested_list_chain, back_translate_nested_list_poly, back_translate_nested_list, back_translate_nested_list_poly_binom,back_translate_nested_list_real,back_translate_nested_list_chain_real
from .oligopool import OligoPool
from .auxiliary import flatten_at_layer


class Assembly(SimulateSynthesis):
    """
    Simulate the assembly of single stranded oligos into long chaines of double stranded DNA.
    """
    def __init__(self, params, logger=None):
        self.params = params
    

    def simulate(self, data):
        """
        Simulate errors that occur when hybridizing and ligation oligos into long chains of double stranded DNA.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        try: 
            synthesis_module = importlib.import_module(f'dnabyte.encoding.{self.params.encoding_method}.assembly')
            synthezised, info = synthesis_module.assembly(data, self.params)

        except KeyError:
            raise ValueError(f"Synthesis method '{self.params.encoding_method}' not recognized.")
        
        return synthezised, info

def attributes(params):
    if 'library' not in params.__dict__ or params.library is None:
        raise ValueError("Library must be specified for synthesis assembly.")
    else:
        library = params.library

    if 'encoding_method' not in params.__dict__ or params.encoding_method is None:
        raise ValueError("Encoding scheme must be specified for synthesis assembly.")
    else:
        encoding_method = params.encoding_method

    if 'assembly_structure' not in params.__dict__ or params.assembly_structure is None:
        raise ValueError("Assembly structure must be specified for synthesis assembly.")
    else:
        assembly_structure = params.assembly_structure

    if 'theory' not in params.__dict__ or params.theory is None:
        theory = 'no'
    else:
        theory = params.theory

    if 'mean' not in params.__dict__ or params.mean is None:
        mean = 10
    else:
        mean = params.mean

    if 'std_dev' not in params.__dict__ or params.std_dev is None:
        std_dev = 0
    else:
        std_dev = params.std_dev
    
    return {
        "library": library,
        "encoding_method": encoding_method,
        "assembly_structure": assembly_structure,
        "theory": theory,
        "mean": mean,
        "std_dev": std_dev
    }