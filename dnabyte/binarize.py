from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode

class Binarize:
    """
    This class handles the binarization and debinarization of data.
    Attributes:
        params: Configuration parameters for binarization.
        name: The name of the binarization method to use.
        binarize_plugins: A dictionary mapping binarization method names to their corresponding classes.
    Methods:
        binarize(data): Binarizes the given data using the specified method.
        debinarize(data, output_file_path=None): Debinarizes the given binary data using the specified method.  
    """
    def __init__(self, params):
        self.params = params
        self.name = params.binarization_method
        self.binarize_plugins = params.binarization_plugins

        if hasattr(params, 'text_encoding'):
            self.text_encoding = params.text_encoding

    def binarize(self, data):
        if self.name == None:
            raise ValueError("No binarization method specified.")
        elif isinstance(data, Data) == False:
            raise TypeError(f"Expected Data object, got {type(data).__name__}")
        else:
            try:
                print(self.binarize_plugins, 'binarize plugins in binarization method')
                binarization_class = self.binarize_plugins[self.name]
                plugin = binarization_class(self.params)  # Instantiate the plugin class

                # Call the binarize method of the plugin class
                print('binarize method called in binarization class')
                successfull, size = plugin.binarize(data)
                print(successfull,data.file_paths, size, 'data in binarization method')
                binary_data = BinaryCode(data=successfull, 
                                     file_paths=data.file_paths,size=size)
                return binary_data
            
                
            except KeyError:
                raise ValueError(f"Binarization method '{self.name}' not recognized.")
            
    def debinarize(self, data, output_file_path=None):
        if self.name == None:
            raise ValueError("No binarization method specified.")
        elif isinstance(data, BinaryCode) == False:
            raise TypeError(f"Expected BinaryCode object, got {type(data).__name__}")
        else:
            try:
                binarization_class = self.binarize_plugins[self.name]
                plugin = binarization_class(self.params)  # Instantiate the plugin class

                # Call the debinarize method of the plugin class
                successful = plugin.debinarize(data, output_file_path)
                return successful

            except KeyError:
                raise ValueError(f"Binarization method '{self.name}' not recognized.")