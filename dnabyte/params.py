from dnabyte.library import Library
from typing import List, Tuple
from pyxdameraulevenshtein import damerau_levenshtein_distance


class Params:
    def __init__(self, 
                 name=None, 
                 filename=None, 
                 assembly_structure=None, 
                 encoding_scheme=None, 
                 inner_error_correction=None, 
                 outer_error_correction=None,
                 dna_barcode_length=None, 
                 codeword_maxlength_positions=None, 
                 sigma_amount=None, 
                 codeword_length=None, 
                 percent_of_symbols=None,
                 index_carry_length=None, 
                 reed_solo_percentage=None, 
                 synthesis_method=None, 
                 mean=None, 
                 vol=None, 
                 std_dev=None, 
                 hybridisation_steps=None, 
                 years=None, 
                 storage_conditions=None, 
                 sequencing_method=None, 
                 iid_error_rate=None,
                 library_name=None,
                 seed=None, 
                 debug=False,
                 theory=None,
                 ltcode_header=None):
        self.name = name
        self.filename = filename
        self.assembly_structure = assembly_structure
        self.encoding_scheme = encoding_scheme
        self.inner_error_correction = inner_error_correction
        self.outer_error_correction = outer_error_correction
        self.dna_barcode_length = dna_barcode_length
        self.codeword_maxlength_positions = codeword_maxlength_positions
        self.sigma_amount = sigma_amount
        self.codeword_length = codeword_length
        self.percent_of_symbols = percent_of_symbols
        self.index_carry_length = index_carry_length
        self.reed_solo_percentage = reed_solo_percentage
        self.synthesis_method = synthesis_method
        self.mean = mean
        self.vol = vol
        self.std_dev = std_dev
        self.hybridisation_steps = hybridisation_steps
        self.years = years
        self.storage_conditions = storage_conditions
        self.sequencing_method = sequencing_method
        self.library_name = library_name
        self.seed = seed
        self.debug = debug
        self.theory = theory
        self.ltcode_header = ltcode_header
        
        # check the validity of the parameters
        if self.sequencing_method == 'iid':
            self.iid_error_rate = iid_error_rate
        
        # check that the library is compatible with the selected assembly structure
        if self.library_name is not None:
            with open('./app/static/libraries/' + self.library_name, 'r') as f:
                first_line = f.readline().strip() 
            if assembly_structure == 'linear_assembly':
                if first_line == 'Messages':
                    raise ValueError("Library not compatible with linear assembly")
                else:
                    self.library = Library(structure=assembly_structure, filename='./app/static/libraries/' + library_name)
            elif assembly_structure == 'positional_assembly':
                if first_line == 'Messages':
                    self.library = Library(structure=assembly_structure, filename='./app/static/libraries/' + library_name)
                else:
                    raise ValueError("Library not compatible with positional assembly")
            else:
                raise ValueError("Library not compatible with synthesis")
            
        #assure message carring part is at least 1 oligo long
        if assembly_structure == 'linear_assembly':
            if encoding_scheme == 'linear_encoding':
                max_codeword_length = codeword_length
            elif encoding_scheme == 'binomial_encoding':
                max_codeword_length = len(self.library.leftmotives) + len(self.library.rightmotives) - 1
        elif assembly_structure == 'positional_assembly':
            max_codeword_length = len(self.library.position) - 1
        if inner_error_correction == 'ltcode':
            message_length = max_codeword_length - dna_barcode_length - codeword_maxlength_positions - 1 - index_carry_length - ltcode_header
        else:
            message_length = max_codeword_length - dna_barcode_length - codeword_maxlength_positions - 1
        if message_length < 1:
            raise ValueError("Invalid code word split.")
                
        if self.encoding_scheme not in ['linear_encoding', 'binomial_encoding', 'max_density_encoding', 'no_homopolymer_encoding','no_homopolymeroddeven_encoding']:
            raise ValueError("Invalid encoding scheme")
        
        if self.assembly_structure not in ['synthesis', 'linear_assembly', 'positional_assembly']:
            raise ValueError("Invalid library structure")
        

        if self.inner_error_correction not in ['ltcode', None]:
            raise ValueError("Invalid inner error correction")
        
        if self.outer_error_correction not in ['reedsolomon', None]:
            raise ValueError("Invalid outer error correction")

        if self.outer_error_correction == 'reedsolomon':
            if self.reed_solo_percentage < 0 or self.reed_solo_percentage > 100:
                raise ValueError("Reedsolomon percentage not compatible")
            
        if self.inner_error_correction == 'ltcode':
            if self.percent_of_symbols < 0 or self.percent_of_symbols > 100:
                raise ValueError("Percent of symbols not compatible")
            
            if self.index_carry_length < 0:
                raise ValueError("Index carry length not compatible")    
        
        if self.assembly_structure == 'positional_assembly' or self.assembly_structure == 'linear_assembly':
            if self.library is None or self.library_name is None:
                raise ValueError("Library not provided")
            
        if self.encoding_scheme == 'binomial_encoding':
            if self.sigma_amount is None or self.theory is None:
                raise ValueError("Missing parameters for encoding scheme")
            
        if self.encoding_scheme == 'linear_encoding' and self.assembly_structure == 'linear_assembly':
            if self.codeword_length is None or self.theory is None:
                raise ValueError("Missing parameters for encoding scheme")
            
        if (self.assembly_structure == 'synthesis' or (self.assembly_structure == 'linear_assembly' and self.encoding_scheme == 'linear_encoding')):
            if self.codeword_length is None:
                raise ValueError("Missing parameters for encoding scheme")
            
        if self.encoding_scheme == 'linear_encoding' and self.assembly_structure == 'positional_assembly':
            if self.theory is None:
                raise ValueError("Missing parameters for encoding scheme")

        if self.sequencing_method not in [41, 40, 37, 36, 39, 38, 35,'iid', None]:
            raise ValueError("Invalid sequencing method")
        
        if self.synthesis_method not in [3, 4, 5, 6, 7, 68, 69, 70, 71, 'nosynthpoly', None]:
            raise ValueError("Invalid Synthesis Error")
        
        if self.storage_conditions not in ['permafrost', 'room_temperature', 'biogene', None]:
            raise ValueError("Invalid storage conditions")

        # NOTE: This is here for testing library distances, to be removed after simulations are done

        # min_distance = float('inf')  # Initialize with infinity
        # closest_strings = []         # List to store strings with the smallest distance

        # for m in self.library.library:
        #     for s in self.library.library:
        #         if m != s:  # Avoid comparing the same string
        #             distance = damerau_levenshtein_distance(m, s)  # Calculate the distance
        #             if distance < min_distance:
        #                 min_distance = distance  # Update the smallest distance
                    
        # print(f"Minimum distance: {min_distance}")
        # breakpoint()

        
    def __str__(self):
        output = f"Name: {self.name}\n"
        output += f"Filename: {self.filename}\n"
        output += f"Assembly structure: {self.assembly_structure}\n"
        output += f"Encoding scheme: {self.encoding_scheme}\n"
        output += f"Inner error correction: {self.inner_error_correction}\n"
        output += f"Outer error correction: {self.outer_error_correction}\n"
        output += f"DNA barcode length: {self.dna_barcode_length}\n"
        output += f"Codeword max length positions: {self.codeword_maxlength_positions}\n"
        output += f"Sigma amount: {self.sigma_amount}\n"
        output += f"Codeword length: {self.codeword_length}\n"
        output += f"Percent of symbols: {self.percent_of_symbols}\n"
        output += f"Index carry length: {self.index_carry_length}\n"
        output += f"Reed-Solo percentage: {self.reed_solo_percentage}\n"
        output += f"Synthesis method: {self.synthesis_method}\n"
        output += f"Mean: {self.mean}\n"
        output += f"Volume: {self.vol}\n"
        output += f"Standard deviation: {self.std_dev}\n"
        output += f"Hybridisation steps: {self.hybridisation_steps}\n"
        output += f"Years: {self.years}\n"
        output += f"Storage conditions: {self.storage_conditions}\n"
        output += f"Sequencing method: {self.sequencing_method}\n"
        output += f"Library name: {self.library_name}\n"
        output += f"Seed: {self.seed}\n"
        output += f"Debug: {self.debug}\n"
        output += f"Theory: {self.theory}\n"
        return output


    @classmethod
    def params_range(cls, **kwargs):
        """
        Alternative constructor that creates a list of Params instances for each item in the list if one of the parameters is a non-empty list.
        
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



        # if params.encoding_scheme in ['linear_encoding', 'binomial_encoding', 'max_density_encoding', 'no_homopolymer_encoding','no_homopolymeroddeven_encoding']:

        #     self.encoding_scheme = params.encoding_scheme
        # else:
        #     raise ValueError("Invalid encoding scheme")
        
        # if params.assembly_structure in ['synthesis', 'linear_assembly', 'positional_assembly']:
        #     self.assembly_structure = params.assembly_structure

        # else:
        #     raise ValueError("Invalid library structure")
            
        # # TODO: check whether the library is compatible with the selected encoding scheme and structure
        # if params.outer_error_correction == 'reedsolomon':
            
        #     if 'reed_solo_percentage' in kwargs:
        #         self.reed_solo_percentage = kwargs['reed_solo_percentage']
        #     else:
        #         raise ValueError("Reedsolomon percentage not provided")
            
        # if assembly_structure == 'positional_assembly' or assembly_structure == 'linear_assembly':
            
        #     self.library = kwargs['library']
        #     self.library_name = kwargs['library_name']
        # else:
        #     self.library = assembly_structure   
        
        # self.inner_error_correction = inner_error_correction
        # self.outer_error_correction = outer_error_correction
        # self.dna_barcode_length = dna_barcode_length
        # self.codeword_maxlength_positions = codeword_maxlength_positions

        # if encoding_scheme == 'binomial_encoding':
        #     if 'sigma_amount' in kwargs:
        #         self.sigma_amount = kwargs['sigmaamount']
        #     else: 
        #         ValueError("Missing parameters for encoding scheme")

        # if 'theory' in kwargs:
        #     self.theory = kwargs['theory']

        # # TODO: check whether all necessary parameters are provided for each given set of encoding scheme and structure

        # if encoding_scheme == 'binomial_encoding' and 'sigmaamount' in kwargs and 'theory' in kwargs:
        #     self.sigma_amount = kwargs['sigma_amount']
        #     self.theory = kwargs['theory']
        # else: 
        #     ValueError("Missing parameters for encoding scheme")

        # if encoding_scheme == 'linear_encoding' and assembly_structure == 'linear_assembly' and 'codeword_length' in kwargs and 'theory' in kwargs:
        #     self.codeword_length = kwargs['codeword_length']
        #     self.theory = kwargs['theory']
        # else:
        #     ValueError("Missing parameters for encoding scheme")

        # #if encoding_scheme == 'binomial_encoding':

        # if (assembly_structure == 'synthesis' or (assembly_structure == 'linear_assembly' and encoding_scheme == 'linear_encoding')) and 'codeword_length' in kwargs:
        #     self.codeword_length = kwargs['codeword_length']
        # else:
        #     ValueError("Missing parameters for encoding scheme")
        
        # if encoding_scheme == 'linear_encoding' and assembly_structure == 'positional_assembly' and 'theory' in kwargs:
        #     self.theory = kwargs['theory']
        # else:
        #     ValueError("Missing parameters for encoding scheme")

        # if inner_error_correction == 'ltcode' and 'percent_of_symbols' in kwargs and 'index_carry_length' in kwargs:
        #     self.percent_of_symbols = kwargs['percent_of_symbols']
        #     self.index_carry_length = kwargs['index_carry_length']
        # else:
        #     ValueError("Missing parameters for inner error correction")
