
from dnabyte.synthesize import SimulateSynthesis
from .oligo import Oligo, complement, translate_nested_list, translate_nested_list_poly,back_translate_nested_list_poly_binom_real, back_translate_nested_list_chain, back_translate_nested_list_poly, back_translate_nested_list, back_translate_nested_list_poly_binom,back_translate_nested_list_real,back_translate_nested_list_chain_real
from .oligopool import OligoPool
from .auxiliary import flatten_at_layer


class Assembly(SimulateSynthesis):
    """
    Simulate the assembly of single stranded oligos into long chaines of double stranded DNA.
    """

    def simulate(self, data):
        """
        Simulate errors that occur when hybridizing and ligation oligos into long chains of double stranded DNA.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """

        # linear assembly / linear encoding
        if self.assembly_structure == 'linear_assembly' and  self.encoding_scheme == 'linear_encoding':
            
            restructuredlib = translate_nested_list(data, self.library.translationlibleft, self.library.translationlibright)
            
            if self.theory == 'yes':
                poolingfinish = self.call_function_repeatedly_therory(restructuredlib, self.library, 'linear_encoding')
                retranslated = back_translate_nested_list_chain(poolingfinish, self.library.translationlibleft, self.library.translationlibright)

            if self.theory == 'no':
                empyrestructuredlib = []
                for i in range(len(restructuredlib)):
                    empyrestructuredlib.append(self.nest_list(self.add_empty_to_innermost(restructuredlib[i], self.library)))
                    
                poolingfinish = self.call_function_repeatedly(empyrestructuredlib, self.library, 'linear_encoding')
                poolingfinish = flatten_at_layer(poolingfinish, 1)
                poolingfinish = OligoPool([]).join(pools=poolingfinish, mean=1, std_dev=0)
                retranslated = back_translate_nested_list_chain_real(poolingfinish, self.library.translationlibleft, self.library.translationlibright)

        # linear assembly / binomial encoding
        elif self.assembly_structure == 'linear_assembly' and  self.encoding_scheme == 'binomial_encoding':

            restructuredlib = translate_nested_list(data, self.library.translationlibleft, self.library.translationlibright)
            
            if self.theory == 'yes':
                poolingfinish = self.call_function_repeatedly_therory(restructuredlib, self.library, 'binomial_encoding')
                retranslated = back_translate_nested_list(poolingfinish, self.library.translationlibleft, self.library.translationlibright)
                retranslated = flatten_at_layer(retranslated, 1)

            elif self.theory == 'no':
                poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'binomial_encoding')
                retranslated = back_translate_nested_list_real(poolingfinish, self.library.translationlibleft, self.library.translationlibright)
                
        # positional assembly / linear encoding
        elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'linear_encoding':
            restructuredlib = translate_nested_list_poly(data, self.library)
            
            if self.theory == 'yes':
                poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'linear_encoding')
                retranslated = back_translate_nested_list_poly(poolingfinish, self.library)
            if self.theory == 'no':
                poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'linear_encoding')
                retranslated = back_translate_nested_list_poly_binom_real(poolingfinish, self.library)

        # positional assembly / binomial encoding
        elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'binomial_encoding':
            restructuredlib = translate_nested_list_poly(data, self.library)
            if self.theory == 'yes':
                poolingfinish = self.call_function_repeatedly_therory(restructuredlib, self.library, 'binomial_encoding')
                retranslated = back_translate_nested_list_poly_binom(poolingfinish, self.library)
                retranslated = flatten_at_layer(retranslated, 1)

            elif self.theory == 'no':
                poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'binomial_encoding')
                retranslated = back_translate_nested_list_poly_binom_real(poolingfinish, self.library)
        
        info = {}

        return retranslated, info

def attributes(params):
    if 'library' not in params.__dict__ or params.library is None:
        raise ValueError("Library must be specified for synthesis assembly.")
    else:
        library = params.library

    if 'encoding_scheme' not in params.__dict__ or params.encoding_scheme is None:
        raise ValueError("Encoding scheme must be specified for synthesis assembly.")
    else:
        encoding_scheme = params.encoding_scheme

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
        "encoding_scheme": encoding_scheme,
        "assembly_structure": assembly_structure,
        "theory": theory
    }