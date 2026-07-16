"""
Post-sequencing processing module using plugin-based architecture.

Routes sequencing output through appropriate processing pipeline
based on encoding method and available post_sequencing_plugins.
"""

from typing import Tuple, Dict, Any
from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.data_classes.clusteredcode import ClusteredDNA


class Cluster:

    
    def __init__(self, params, logger=None):
        self.clustering_method = params.clustering_method
        self.clustering_plugins = params.clustering_plugins
        self.params = params
        self.logger = logger
        
    def cluster(self, data):
        """
        Simulate clustering of DNA sequences.
        
        :param data: An object of class InSilicoDNA.
        :return: A list of clustered DNA sequences.
        """
        
        if isinstance(data, InSilicoDNA):
            # Dynamically find the appropriate class based on clustering_method
            try:
                data_to_cluster = data.data
                clustering_class = self.clustering_plugins[self.clustering_method.lower()]
                for key, value in self.clustering_plugins.items():
                    setattr(self, key, value)
                plugin = clustering_class(self.params, self.logger)  # Instantiate the plugin class
                data_sto, info = plugin.cluster(data_to_cluster)
                obj = ClusteredDNA(data_sto)
                if hasattr(data, 'file_paths'):
                    obj.file_paths = data.file_paths  
                return obj, info
            except KeyError:
                raise ValueError(f"Clustering method '{self.clustering_method}' is not recognized. ")
        else:
            raise ValueError("The input data is not an instance of InSilicoDNA.")        