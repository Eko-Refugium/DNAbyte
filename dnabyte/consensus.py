
from typing import Tuple, Dict, Any
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.data_classes.clusteredcode import ClusteredDNA


class Consensus:

    
    def __init__(self, params, logger=None):
        self.consensus_method = params.recovery_method
        self.consensus_plugins = params.recovery_plugins
        self.params = params
        self.logger = logger
        
    def call(self, data):
        """
        Simulate consensus of DNA sequences.
        
        :param data: An object of class InSilicoDNA.
        :return: A list of consensus DNA sequences.
        """
        
        if isinstance(data, ClusteredDNA):
            # Dynamically find the appropriate class based on consensus_method
            try:
                data_to_consensus = data.data
                print(f"Data to consensus: {self.consensus_plugins}")
                consensus_class = self.consensus_plugins[self.consensus_method.lower()]
                for key, value in self.consensus_plugins.items():
                    setattr(self, key, value)
                print(f"Initializing consensus plugin for method: {self.consensus_method}")
                plugin = consensus_class(self.params, self.logger)  # Instantiate the plugin class
                print(f"Using consensus plugin: {self.consensus_method}")
                data_sto, info = plugin.recover(data_to_consensus)
                print("RECOVER RETURN TYPE:", type(data_sto))
                print("RECOVER RETURN VALUE:", data_sto)
                obj = NucleobaseCode(data_sto)
                if hasattr(data, 'file_paths'):
                    obj.file_paths = data.file_paths
                return obj, info
            except KeyError:
                raise ValueError(f"Consensus method '{self.consensus_method}' is not recognized. ")
        else:
            raise ValueError("The input data is not an instance of ClusteredDNA.")       