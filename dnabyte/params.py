import importlib
from dnabyte.load_plugins import load_plugins

class Params:

    def __init__(self, debug=False, **kwargs):

        # Set default values for required attributes
        self.encoding_method = None
        self.synthesis_method = None
        self.storage_conditions = None
        self.sequencing_method = None
        self.binarization_method = None
        
        # Set parameters from kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.debug = debug

        print(self)

        # Load plugins
        self.binarization_plugins, self.encoding_plugins, self.synthesis_plugins, self.storage_plugins, self.sequencing_plugins  = load_plugins(self.encoding_method, self.synthesis_method, self.storage_conditions, self.sequencing_method, self.binarization_method)

        # Check binarization parameters
        if self.binarization_method is None:
            pass
        elif self.binarization_method in self.binarization_plugins:
            binarization = importlib.import_module(f"dnabyte.binarization.{self.binarization_method}.binarize")
            attributes_bin = binarization.attributes(self)
            for keys, value in attributes_bin.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid binarization method: {self.binarization_method}")
        
        # Check encoding parameters
        if self.encoding_method is None:
            pass
        elif self.encoding_method in self.encoding_plugins:
            encoding = importlib.import_module(f"dnabyte.encoding.{self.encoding_method}.encode")
            attributes_enc = encoding.attributes(self) 
            for keys, value in attributes_enc.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid encoding method: {self.encoding_method}")
        
        # Check synthesis parameters
        if self.synthesis_method is None:
            pass
        elif 'synthesis_' + self.synthesis_method in self.synthesis_plugins:
            synthesis = importlib.import_module(f"dnabyte.synthesis.{self.synthesis_method}.synthesize")
            print("Synthesis module loaded:", synthesis)
            attributes_synth = synthesis.attributes(self)
            print("Synthesis attributes:", attributes_synth)
            for keys, value in attributes_synth.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid synthesis method: {self.synthesis_method}")
        
        # Check storage parameters
        if self.storage_conditions is None:
            pass
        elif 'storage_' + self.storage_conditions in self.storage_plugins:
            storage = importlib.import_module(f"dnabyte.storage.{self.storage_conditions}.store")
            attributes_store = storage.attributes(self)
            for keys, value in attributes_store.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid storage conditions: {self.storage_conditions}")
        
        # Check sequencing parameters
        if self.sequencing_method is None:
            pass
        elif 'sequencing_' + self.sequencing_method in self.sequencing_plugins:
            sequencing = importlib.import_module(f"dnabyte.sequencing.{self.sequencing_method}.sequence")
            attributes_seq = sequencing.attributes(self)
            for keys, value in attributes_seq.items():
                setattr(self, keys, value)
        else:
            raise ValueError(f"Invalid sequencing method: {self.sequencing_method}")
        


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

