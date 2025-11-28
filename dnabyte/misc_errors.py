import os
import importlib
import inspect
import random
import math

from dnabyte.data_classes.insilicodna import InSilicoDNA

class SimulateMiscErrors:
    
    def __init__(self, params, logger=None):
        self.error_para = params.error_para
        self.error_conditions = params.error_conditions
        self.error_plugins = params.error_plugins

    def simulate(self, to_error_data):
        
        if self.error_conditions == None:
            stored_data = InSilicoDNA(to_error_data.data)
            info = {"number_of_strand_breaks": 0}
            return stored_data, info
        
        elif isinstance(to_error_data, InSilicoDNA):
            # Dynamically find the appropriate class based on storage_conditions
            try:
                if isinstance(self.error_conditions, str):
                    print(self.error_plugins, 'error conditions in simulate misc errors')
                    error_class = self.error_plugins['error_' + self.error_conditions.lower()]
                    print(error_class, 'error class in simulate misc errors')
                    plugin = error_class(self)  # Instantiate the plugin class
                    data_sto, info = plugin.simulate(to_error_data.data)
                    data_sto = InSilicoDNA(data_sto)
                    return data_sto, info
                elif isinstance(self.error_conditions, list):
                    combined_data = to_error_data.data
                    combined_info = {}
                    tempconditions = self.error_para
                    for i, condition in enumerate(self.error_conditions):
                        print(self.error_plugins, 'error conditions in simulate misc errors - list')
                        error_class = self.error_plugins['error_' + condition.lower()]
                        print(error_class, 'error class in simulate misc errors - list')
                        self.error_para = tempconditions[i]
                        plugin = error_class(self)  # Instantiate the plugin class
                        combined_data, info = plugin.simulate(combined_data)
                    return InSilicoDNA(combined_data), info
            except KeyError:
                raise ValueError(f"Error condition '{self.error_conditions}' is not recognized. ")
        else:
            raise ValueError("The input data is not an instance of InSilicoDNA.")        
