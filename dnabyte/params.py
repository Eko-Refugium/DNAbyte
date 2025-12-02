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
    def __init__(self, 
                 name=None, 
                 filename=None, 
                 assembly_structure=None, 
                 encoding_scheme=None, 
                 encoding_method=None,
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
                 mesa_id=3,
                 mean=None, 
                 vol=None, 
                 std_dev=None, 
                 hybridisation_steps=None, 
                 years=None, 
                 storage_conditions=None, 
                 sequencing_method=None, 
                 mesa_sequencing_id=None,
                 iid_error_rate=None,
                 library_name=None,
                 seed=None, 
                 debug=False,
                 theory='no',
                 ltcode_header=None,
                 binarization_method=None,
                 text_encoding='utf-8'):
        self.name = name
        self.filename = filename
        self.assembly_structure = assembly_structure
        self.encoding_scheme = encoding_scheme
        self.encoding_method = encoding_method
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
        self.mesa_id = mesa_id
        self.mean = mean
        self.vol = vol
        self.std_dev = std_dev
        self.hybridisation_steps = hybridisation_steps
        self.years = years
        self.storage_conditions = storage_conditions
        self.sequencing_method = sequencing_method
        self.mesa_sequencing_id = mesa_sequencing_id
        self.library_name = library_name
        self.seed = seed
        self.debug = debug
        self.theory = theory
        self.ltcode_header = ltcode_header
        self.binarization_method = binarization_method
        self.text_encoding = text_encoding
        
        # Load plugins
        self.encoding_plugins, self.synthesis_plugins, self.storage_plugins, self.sequencing_plugins, self.binarization_plugins = load_plugins(self.encoding_method, self.synthesis_method, self.storage_conditions, self.sequencing_method, self.binarization_method)
        # print(self.encoding_plugins)
        # print(self.synthesis_plugins)
        # print(self.storage_plugins)
        # print(self.sequencing_plugins)
        # print(self.binarization_plugins)
        
        if self.mean is None:
            self.mean = 10
        if self.std_dev is None:
            self.std_dev = 2
        if self.vol is None:
            self.vol = 1e-6
        if self.hybridisation_steps is None:
            self.hybridisation_steps = 1000
        if self.years is None:
            self.years = 0

        print(self.encoding_method,self.encoding_plugins, 'encoding method in params')
        
        if self.encoding_method is None:
            self.assembly_structure = None
            self.library = None
            self.encoding_structure = None
        elif self.encoding_method in self.encoding_plugins:
            encoding = importlib.import_module(f"dnabyte.encoding.{self.encoding_method}.encode")
            self.library, self.library_name, self.encoding_scheme, self.assembly_structure, self.dna_barcode_length, self.codeword_maxlength_positions, self.codeword_length, self.percent_of_symbols, self.index_carry_length, self.sigma_amount, self.ltcode_header, self.reed_solo_percentage = encoding.attributes(self)
        else:
            raise ValueError(f"Invalid encoding method: {self.encoding_method}")

        # check the validity of the parameters
        if self.sequencing_method == 'iid':
            self.iid_error_rate = iid_error_rate
        
                
        # if self.encoding_method not in ['linear_binom', 'linear_chain', 'max_density', 'nohomopolymer','polybinom', 'polychain', None]:
        #     raise ValueError("Invalid encoding scheme")
        
        #with new library structure NOT NEEDED
        # if self.assembly_structure not in ['synthesis', 'linear_assembly', 'positional_assembly']:
        #     raise ValueError("Invalid library structure")
        

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

        if self.sequencing_method not in ['mesa', 'iid', 'illumina', 'nanopore', None]:
            raise ValueError("Invalid sequencing method")
        
        if self.mesa_sequencing_id not in [41, 40, 37, 36, 39, 38, 35,'iid', None]:
            raise ValueError("Invalid sequencing method")
        
        if self.synthesis_method not in ['assembly', 'mesa', 'nosynthpoly', None]:
            raise ValueError("Invalid synthesis method")
        
        if self.mesa_id not in [3, 4, 5, 6, 7, 68, 69, 70, 71]:
            raise ValueError("Invalid MESA ID")
        
        # Validate storage conditions
        if isinstance(self.storage_conditions, list):
            for condition in self.storage_conditions:
                if condition not in self.storage_plugins and condition not in ['permafrost', 'roomtemperature', 'biogene', None]:
                    raise ValueError("Invalid storage conditions")
        elif self.storage_conditions is None:
            pass
        elif isinstance(self.storage_conditions, str):
            if self.storage_conditions not in self.storage_plugins and self.storage_conditions not in ['permafrost', 'roomtemperature', 'biogene', None]:
                raise ValueError("Invalid storage conditions")
        else:
            raise ValueError("Invalid storage conditions")

        # if self.binarization_method == None:
        #     self.binarization_method = 'binarize_default'
        # print(self.binarization_method, self.binarization_plugins, 'binarization method in params')
        # if self.binarization_method not in self.binarization_plugins:
        #     raise ValueError("Invalid binarization method")


        # # Validate storage conditions
        # if self.storage_conditions and self.storage_conditions.lower() not in self.storage_plugins:
        #     raise ValueError(
        #         f"Invalid storage condition '{self.storage_conditions}'. "
        #         f"Available options are: {', '.join(self.storage_plugins.keys())}."
        #     )

        # # Validate sequencing methods
        # if self.sequencing_method and self.sequencing_method.lower() not in self.sequencing_plugins:
        #     raise ValueError(
        #         f"Invalid sequencing method '{self.sequencing_method}'. "
        #         f"Available options are: {', '.join(self.sequencing_plugins.keys())}."
        #     )


    # def load_plugins(self):
    #     """
    #     Load plugins from the specified folder and return a dictionary of plugin classes.
    #     """
    #     plugin_folder = os.path.dirname(__file__)
    #     encoding_plugins = {}  # Define the plugins dictionary
    #     synthesis_plugins = {}  # Define the plugins dictionary
    #     storage_plugins = {}  # Define the plugins dictionary
    #     sequencing_plugins = {}  # Define the plugins dictionary
    #     binarization_plugins = {}  # Define the plugins dictionary

    #     path_encode = plugin_folder + '/EncodingMethods'

    #     folders_encode = [name for name in os.listdir(path_encode) 
    #         if os.path.isdir(os.path.join(path_encode, name))]
        
    #     # path_sequencing = plugin_folder + '/ErrorChannels/Sequencing'
    #     # folders_sequencing = [name for name in os.listdir(path_sequencing) 
    #     #     if os.path.isdir(os.path.join(path_sequencing, name))]
        
    #     # path_synthesis = plugin_folder + '/ErrorChannels/Synthesis'
    #     # folders_synthesis = [name for name in os.listdir(path_synthesis) 
    #     #     if os.path.isdir(os.path.join(path_synthesis, name))]   
        
    #     # path_sortage = plugin_folder + '/ErrorChannels/Storage'
    #     # folders_storage = [name for name in os.listdir(path_sortage)
    #     #     if os.path.isdir(os.path.join(path_sortage, name))]
        
    #     # path_errorchannels = plugin_folder + '/ErrorCorrection'
    #     # folders_errorcorrection = [name for name in os.listdir(path_errorchannels)
    #     #     if os.path.isdir(os.path.join(path_errorchannels, name))]

    #     # return folders_encode, folders_synthesis, folders_storage, folders_sequencing, folders_errorcorrection

    #     # Load encoding plugins
    #     for filename in folders_encode:
    #         # if filename.startswith('encode_') and filename.endswith('.py'):
    #         # encoding_type = filename[7:-3]  # Extract the encoding type (e.g., 'maxdensity' from 'encode_maxdensity.py')
    #         try:
    #             # Load the encode module
    #             #print('dnabyte.EncodingMethods' + f'.{filename}.encode')
    #             encode_module = importlib.import_module('dnabyte.EncodingMethods' + f'.{filename}.encode')
    #             #print(encode_module)
    #             # Find the class definition in the module
    #             for name, obj in inspect.getmembers(encode_module, inspect.isclass):
    #                 if issubclass(obj, Encode) and obj is not Encode:
    #                     # Add the class to the plugins dictionary
    #                     encoding_plugins[filename] = obj
    #                     break
    #         except Exception as e:
    #             print(f"Error loading encoding plugin '{filename}': {e}")
        
    #     for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/ErrorChannels/Sequencing') if f.endswith('.py')]:
    #         # if filename.startswith('encode_') and filename.endswith('.py'):
    #         module_name = filename[:-3] 
    #         try:
    #             # Load the encode module
    #             sequencing_module = importlib.import_module('dnabyte.ErrorChannels.Sequencing' + f'.{module_name}')
    #             # Find the class definition in the module
    #             for name, obj in inspect.getmembers(sequencing_module, inspect.isclass):
    #                 if issubclass(obj, SimulateSequencing) and obj is not SimulateSequencing:
    #                     # Add the class to the plugins dictionary
    #                     sequencing_plugins[module_name] = obj
    #                     break
    #         except Exception as e:
    #             print(f"Error loading sequencing plugin '{module_name}': {e}")

    #     for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/ErrorChannels/Synthesis') if f.endswith('.py')]:
    #         # if filename.startswith('encode_') and filename.endswith('.py'):
    #         module_name = filename[:-3] 
    #         try:
    #             # Load the encode module
    #             synthesis_module = importlib.import_module('dnabyte.ErrorChannels.Synthesis' + f'.{module_name}')
    #             # Find the class definition in the module
    #             for name, obj in inspect.getmembers(synthesis_module, inspect.isclass):
    #                 if issubclass(obj, SimulateAssembly) and obj is not SimulateAssembly:
    #                     # Add the class to the plugins dictionary
    #                     synthesis_plugins[module_name] = obj
    #                     break
    #         except Exception as e:
    #             print(f"Error loading synthesis plugin '{module_name}': {e}")
        
    #     for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/ErrorChannels/Storage') if f.endswith('.py')]:
    #         # if filename.startswith('encode_') and filename.endswith('.py'):
    #         #print(filename)
    #         module_name = filename[:-3]
    #         try:
    #             # Load the encode module
    #             storage_module = importlib.import_module('dnabyte.ErrorChannels.Storage' + f'.{module_name}')
    #             # Find the class definition in the module
    #             for name, obj in inspect.getmembers(storage_module, inspect.isclass):
    #                 if issubclass(obj, SimulateStorage) and obj is not SimulateStorage:
    #                     # Add the class to the plugins dictionary
    #                     storage_plugins[module_name] = obj
    #                     break
    #         except Exception as e:
    #             print(f"Error loading storage plugin '{module_name}': {e}")

    #     for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/binarization') if f.endswith('.py')]:
    #         #print(filename)
    #         module_name = filename[:-3] 
    #         try:
    #             # Load the encode module
    #             binarization_module = importlib.import_module('dnabyte.binarization' + f'.{module_name}')
    #             # Find the class definition in the module
    #             for name, obj in inspect.getmembers(binarization_module, inspect.isclass):
    #                 if issubclass(obj, Binarize) and obj is not Binarize:
    #                     # Add the class to the plugins dictionary
    #                     binarization_plugins[module_name] = obj
    #                     break
    #         except Exception as e:
    #             print(f"Error loading binarization plugin '{module_name}': {e}")

        # Load error channel plugins
        # for filename in :
        #     if filename.endswith('.py') and (filename.startswith('synthesis_') or filename.startswith('storage_') or filename.startswith('sequencing_')):
        #         module_name = filename[:-3]  # Remove the '.py' extension
        #         module = importlib.import_module(f"dnabyte.{module_name}")
        #         for name, obj in inspect.getmembers(module, inspect.isclass):
        #             # Check if the class is a subclass of SimulateSynthesis
        #             if filename.startswith('synthesis_') and issubclass(obj, SimulateSynthesis) and obj is not SimulateSynthesis:
        #                 synthesis_plugins[name.lower()] = obj
        #             # Check if the class is a subclass of SimulateStorage
        #             elif filename.startswith('storage_') and issubclass(obj, SimulateStorage) and obj is not SimulateStorage:
        #                 storage_plugins[name.lower()] = obj
        #             # Check if the class is a subclass of SimulateSequencing
        #             elif filename.startswith('sequencing_') and issubclass(obj, SimulateSequencing) and obj is not SimulateSequencing:
        #                 sequencing_plugins[name.lower()] = obj

    #     return encoding_plugins, synthesis_plugins, storage_plugins, sequencing_plugins, binarization_plugins


    #     # NOTE: This is here for testing library distances, to be removed after simulations are done

    #     # min_distance = float('inf')  # Initialize with infinity
        # closest_strings = []         # List to store strings with the smallest distance

        # for m in self.library.library:
        #     for s in self.library.library:
        #         if m != s:  # Avoid comparing the same string
        #             distance = damerau_levenshtein_distance(m, s)  # Calculate the distance
        #             if distance < min_distance:
        #                 min_distance = distance  # Update the smallest distance

        
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
