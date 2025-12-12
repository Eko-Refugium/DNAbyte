import os
import importlib
import inspect
import random
import math

from dnabyte.data_classes.insilicodna import InSilicoDNA

class SimulateMiscErrors:
    
    def __init__(self, params, logger=None):
        self.error_plugins = params.error_plugins
        if not hasattr(params, 'error_methods') or params.error_methods is None:
            self.error_methods = None
        else:
            self.error_methods = params.error_methods
            if hasattr(params, 'error_params'):
                self.error_params = params.error_params
            elif hasattr(params, 'error_params_list'):
                self.error_params_list = params.error_params_list
            else:
                raise ValueError("Params object must have error_conditions for misc error plugin.")
    def simulate(self, to_error_data):
        
        if self.error_methods == None:
            stored_data = InSilicoDNA(to_error_data.data)
            info = {"number_of_strand_breaks": 0}
            return stored_data, info
        
        elif isinstance(to_error_data, InSilicoDNA):
            # Dynamically find the appropriate class based on storage_conditions
            try:
                if isinstance(self.error_methods, str):
                    error_class = self.error_plugins[self.error_methods.lower()]
                    plugin = error_class(self)  # Instantiate the plugin class
                    data_sto, info = plugin.simulate(to_error_data.data)
                    data_sto = InSilicoDNA(data_sto)
                    return data_sto, info
                elif isinstance(self.error_methods, list):
                    combined_data = to_error_data.data
                    tempconditions = self.error_params_list
                    for i, condition in enumerate(self.error_methods):
                        error_class = self.error_plugins[condition.lower()]
                        self.error_params = tempconditions[condition]
                        for key, value in self.error_params.items():
                            setattr(self, key, value)
                        plugin = error_class(self)  # Instantiate the plugin class
                        combined_data, info = plugin.simulate(combined_data)
                    return InSilicoDNA(combined_data), info
            except KeyError:
                raise ValueError(f"Error condition '{self.error_methods}' is not recognized. ")
        else:
            raise ValueError("The input data is not an instance of InSilicoDNA.")        
