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

    def simulate(self, data):
        """
        Simulate sequencing of DNA sequences.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        if self.sequencing_method == None:
            sequenced_data = InSilicoDNA(data=data.data)
            info = {"error_dict": {}, "number_of_sequencing_errors": 0}
            return sequenced_data, info

        else:
            if isinstance(data, InSilicoDNA):
            # Dynamically find the appropriate class based on sequencing_method
                try:
                    if self.params.sequencing_method in [41, 40, 37, 36, 39, 38, 35]:
                        sequencing_class = self.sequencing_plugins["mesa"]
                        plugin = sequencing_class(self)
                    else:
                        sequencing_class = self.sequencing_plugins[self.sequencing_method.lower()]
                        plugin = sequencing_class(self)  # Instantiate the plugin class

                    # Call the simulate method of the plugin class
                    data_seq, info = plugin.simulate(data)
                    data_seq = InSilicoDNA(data=data_seq)
                    return data_seq, info
                
                except KeyError:
                    raise ValueError(f"Sequencing method '{self.sequencing_method}' not recognized.")
            else:
                raise ValueError("The input data is not an instance of StoredData.")
            
