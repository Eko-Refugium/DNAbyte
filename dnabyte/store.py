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
        
        if not hasattr(params, 'storage_conditions') or params.storage_conditions is None:
            self.storage_conditions = None
        else:
            self.years = params.years
            self.storage_conditions = params.storage_conditions
            self.storage_plugins = params.storage_plugins
            self.logger = logger

    def simulate(self, data):
        """
        Simulate storage of DNA sequences in a storage medium.
        
        :param data: An object of class InSilicoDNA.
        :return: A list of stored DNA sequences.
        """
        
        if isinstance(data, InSilicoDNA):
            # Dynamically find the appropriate class based on storage_conditions
            try:
                if self.storage_conditions == None:
                    stored_data = InSilicoDNA(data=data.data)
                    info = {"degradation_info": {}, "number_of_strand_breaks": 0}
                    return stored_data, info
                elif isinstance(self.storage_conditions, str):
                    storage_class = self.storage_plugins[self.storage_conditions.lower()]
                    plugin = storage_class(self)  # Instantiate the plugin class
                    data_sto, info = plugin.simulate(data)
                    return InSilicoDNA(data_sto), info
                
                elif isinstance(self.storage_conditions, list):
                    combined_data = data.data
                    combined_info = {}
                    tempconditions = self.years
                    for i, condition in enumerate(self.storage_conditions):
                        storage_class = self.storage_plugins[condition.lower()]
                        self.years = tempconditions[i]
                        plugin = storage_class(self)  # Instantiate the plugin class
                        combined_data, info = plugin.simulate(combined_data)
                    return InSilicoDNA(combined_data), info
            except KeyError:
                raise ValueError(f"Storage condition '{self.storage_conditions}' is not recognized. ")
        else:
            raise ValueError("The input data is not an instance of InSilicoDNA.")        
