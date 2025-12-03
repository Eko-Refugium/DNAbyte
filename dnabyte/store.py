import os
import importlib
import inspect
import random
import math

from dnabyte.data_classes.insilicodna import InSilicoDNA

class SimulateStorage:
    """
    Simulate storage of DNA sequences in a storage medium.

    The class consists of a constructor, which sets the parameters of the storage simulation, 
    and a simulate method, which applies the simulation to an object of the AssembledData class.
    
    :param years: The number of years for which DNA storage is supposed to be simulated. (0 to perform no simulation)
    :param storage_conditions: A string 'permafrost', or 'room_temperature', or 'biogene'. (None to perform no simulation)
    """
    
    def __init__(self, params, logger=None):
        if not hasattr(params, 'years') or params.years is None:
            self.years = 0
        else:
            self.years = params.years
        self.storage_conditions = params.storage_conditions
        self.storage_plugins = params.storage_plugins

    def simulate(self, assembled_data):
        """
        Simulate storage of DNA sequences in a storage medium.
        
        :param assembled_data: An object of class AssembledData.
        :return: A list of stored DNA sequences.
        """
        if self.storage_conditions == None or self.years == 0:
            stored_data = InSilicoDNA(assembled_data.data)
            info = {"number_of_strand_breaks": 0}
            return stored_data, info
        
        elif isinstance(assembled_data, InSilicoDNA):
            # Dynamically find the appropriate class based on storage_conditions
            try:
                if isinstance(self.storage_conditions, str):
                    print(self.storage_plugins, 'storage conditions in simulate storage')
                    storage_class = self.storage_plugins[self.storage_conditions.lower()]
                    print(storage_class, 'storage class in simulate storage')
                    plugin = storage_class(self)  # Instantiate the plugin class
                    data_sto, info = plugin.simulate(assembled_data.data)
                    data_sto = InSilicoDNA(data_sto)
                    return data_sto, info
                elif isinstance(self.storage_conditions, list):
                    combined_data = assembled_data.data
                    combined_info = {}
                    tempconditions = self.years
                    for i, condition in enumerate(self.storage_conditions):
                        print(self.storage_plugins, 'storage conditions in simulate storage - list')
                        storage_class = self.storage_plugins[condition.lower()]
                        print(storage_class, 'storage class in simulate storage - list')
                        self.years = tempconditions[i]
                        print(self.years, 'self storage conditions in simulate storage - list')
                        plugin = storage_class(self)  # Instantiate the plugin class
                        combined_data, info = plugin.simulate(combined_data)
                    return InSilicoDNA(combined_data), info
            except KeyError:
                raise ValueError(f"Storage condition '{self.storage_conditions}' is not recognized. ")
        else:
            raise ValueError("The input data is not an instance of InSilicoDNA.")        
