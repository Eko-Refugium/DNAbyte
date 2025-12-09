import numpy as np
import random

from dnabyte.data_classes.insilicodna import InSilicoDNA

class SimulateSequencing:
    """
    Simulate sequencing of DNA sequences.
    There are two main types of sequencing technologies: Illumina and Nanopore, each with their 
    characteristic error profiles.

    :param sequencing_technology: A string 'illumina', or 'nanopore'. None to perform no sequencing simulation.
    """
    def __init__(self, params, logger=None):
        self.sequencing_method = params.sequencing_method
        self.sequencing_plugins = params.sequencing_plugins
        self.params = params
        self.logger = logger

    def simulate(self, data):
        """
        Simulate sequencing of DNA sequences.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """

        if isinstance(data, InSilicoDNA):
        # Dynamically find the appropriate class based on sequencing_method
            try:
                if self.sequencing_method == None:
                    obj = InSilicoDNA(data=data.data)
                    if hasattr(data, 'file_paths'):
                        obj.file_paths = data.file_paths
                    return obj, {}
                else:
                    sequencing_class = self.sequencing_plugins[self.sequencing_method]
                    plugin = sequencing_class(self.params, logger=self.logger)
                    data, info = plugin.simulate(data.data)
                    # print(data)
                    obj = InSilicoDNA(data=data)
                    if hasattr(data, 'file_paths'):
                        obj.file_paths = data.file_paths
                    return obj, info                
            except KeyError:
                raise ValueError(f"Sequencing method '{self.sequencing_method}' not recognized.")
        else:
            raise ValueError("The input data is not an instance of InSilicoDNA.")
            


            # if self.params.sequencing_method in [41, 40, 37, 36, 39, 38, 35]:
            #     sequencing_class = self.sequencing_plugins["mesa"]
            #     plugin = sequencing_class(self)
            # else:
            #     sequencing_class = self.sequencing_plugins[self.sequencing_method.lower()]
            #     plugin = sequencing_class(self)  # Instantiate the plugin class

            # # Call the simulate method of the plugin class
            # data_seq, info = plugin.simulate(data)
            # data_seq = InSilicoDNA(data=data_seq)
            # return data_seq, info
        
