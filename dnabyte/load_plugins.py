import os
import importlib
import inspect

from dnabyte.consensus import Consensus
from dnabyte.encode import Encode
from dnabyte.store import SimulateStorage

from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.misc_err import SimulateMiscErrors
from dnabyte.cluster import Cluster


def load_plugins(binarization_method, encoding_method, storage_conditions, sequencing_method, error_methods, clustering_method=None, recovery_method=None):
    """
    Load plugins from the specified folder and return a dictionary of plugin classes.
    """
    plugin_folder = os.path.dirname(__file__)
    encoding_plugins = {}  # Define the plugins dictionary
    storage_plugins = {}  # Define the plugins dictionary
    sequencing_plugins = {}  # Define the plugins dictionary
    error_plugins = {}  # Define the plugins dictionary
    binarization_plugins = {}  # Define the plugins dictionary
    clustering_plugins = {}  # Define the plugins dictionary
    recovery_plugins = {}  # Define the plugins dictionary

    # Normalize method names to lowercase
    binarization_method = binarization_method.lower() if binarization_method != None else None
    encoding_method = encoding_method.lower() if encoding_method != None else None
    # if isinstance(synthesis_method, str):
    #     synthesis_method = synthesis_method.lower() if synthesis_method != None else None
    if isinstance(storage_conditions, str):
        storage_conditions = storage_conditions.lower() if storage_conditions != None else None
    elif isinstance(storage_conditions, list):
        storage_conditions = [condition.lower() for condition in storage_conditions]
    print(sequencing_method, type(sequencing_method))
    sequencing_method = sequencing_method.lower() if sequencing_method != None else None
        

    # Load binarization plugins
    if binarization_method != None:
        try:
            folders_binarize = [name for name in os.listdir(plugin_folder + '/binarization') 
                if os.path.isdir(os.path.join(plugin_folder + '/binarization', name))]
        except Exception as e:
            print(f"Error listing binarization plugins: {e}")
            folders_binarize = []

        for filename in folders_binarize:
            try:
                if binarization_method.lower() == filename.lower():
                    
                    # Load the binarization module
                    binarization_module = importlib.import_module(f'dnabyte.binarization.{filename}.binarize')
                    # Find the class definition in the module
                    for name, obj in inspect.getmembers(binarization_module, inspect.isclass):
                        if issubclass(obj, Binarize) and obj is not Binarize:
                        # Add the class to the plugins dictionary
                            binarization_plugins[filename] = obj
                            break
            except Exception as e:
                print(f"Error loading binarization plugin '{filename}': {e}")
    

    # Load encoding plugins
    if encoding_method != None:
        try:
            folders_encode = [name for name in os.listdir(plugin_folder + '/encoding') 
                if os.path.isdir(os.path.join(plugin_folder + '/encoding', name))]
        except Exception as e:
            print(f"Error listing encoding plugins: {e}")
            folders_encode = []
        
        for filename in folders_encode:
            try:
                if encoding_method.lower() == filename.lower():
                    # Load the encode module
                    encode_module = importlib.import_module(f'dnabyte.encoding.{filename}.encode')

                    for name, obj in inspect.getmembers(encode_module, inspect.isclass):
                        if issubclass(obj, Encode) and obj is not Encode:
                            # Add the class to the plugins dictionary
                            encoding_plugins[filename] = obj
                            break
            except Exception as e:
                print(f"Error loading encoding plugin '{filename}': {e}")

    # Load storage plugins
    if storage_conditions != None:
        try:   
            folders_store = [name for name in os.listdir(plugin_folder + '/storage') 
                if os.path.isdir(os.path.join(plugin_folder + '/storage', name))]
            
            for filename in folders_store:

                if isinstance(storage_conditions, str):
                    if storage_conditions.lower() == filename.lower():
                        # Load the encode module
                        storage_module = importlib.import_module(f'dnabyte.storage.{filename}.store')
                        # Find the class definition in the module
                        for name, obj in inspect.getmembers(storage_module, inspect.isclass):
                            if issubclass(obj, SimulateStorage) and obj is not SimulateStorage:
                                # Add the class to the plugins dictionary
                                storage_plugins[filename] = obj
                                break
                elif isinstance(storage_conditions, list):
                    for condition in storage_conditions:
                        if condition.lower() == filename.lower():
                            # Load the encode module
                            storage_module = importlib.import_module(f'dnabyte.storage.{filename}.store')
                            # Find the class definition in the module
                            for name, obj in inspect.getmembers(storage_module, inspect.isclass):
                                if issubclass(obj, SimulateStorage) and obj is not SimulateStorage:
                                    # Add the class to the plugins dictionary
                                    storage_plugins[filename] = obj
                                    break
        except Exception as e:
            print(f"Error loading storage plugin '{filename}': {e}")

    # Load error plugins
    if error_methods != None:
        try:   
            folders_error = [name for name in os.listdir(plugin_folder + '/misc_errors') 
                if os.path.isdir(os.path.join(plugin_folder + '/misc_errors', name))]
            for filename in folders_error:

                if isinstance(error_methods, str):
                    if error_methods.lower() == filename.lower():
                        # Load the encode module
                        error_module = importlib.import_module(f'dnabyte.misc_errors.{filename}.err')
                        # Find the class definition in the module
                        for name, obj in inspect.getmembers(error_module, inspect.isclass):
                            if issubclass(obj, SimulateMiscErrors) and obj is not SimulateMiscErrors:
                                # Add the class to the plugins dictionary
                                error_plugins[filename] = obj
                                break
                elif isinstance(error_methods, list):
                    for condition in error_methods:
                        if condition.lower() == filename.lower():
                            # Load the encode module
                            error_module = importlib.import_module(f'dnabyte.misc_errors.{filename}.err')
                            # Find the class definition in the module
                            for name, obj in inspect.getmembers(error_module, inspect.isclass):
                                if issubclass(obj, SimulateMiscErrors) and obj is not SimulateMiscErrors:
                                    # Add the class to the plugins dictionary
                                    error_plugins[filename] = obj
                                    break
        except Exception as e:
            
            print(f"Error loading storage plugin '{filename}': {e}")

    # Load sequencing plugins
    if sequencing_method != None:
        try:
            folders_sequence = [name for name in os.listdir(plugin_folder + '/sequencing') 
                if os.path.isdir(os.path.join(plugin_folder + '/sequencing', name))]
            
            for filename in folders_sequence:
                if sequencing_method.lower() == filename.lower():
                    # Load the encode module
                    sequencing_module = importlib.import_module(f'dnabyte.sequencing.{filename}.sequence')

                    # Find the class definition in the module
                    for name, obj in inspect.getmembers(sequencing_module, inspect.isclass):
                        if issubclass(obj, SimulateSequencing) and obj is not SimulateSequencing:
                            # Add the class to the plugins dictionary
                            sequencing_plugins[filename] = obj
                            break
        except Exception as e:
            print(f"Error loading sequencing plugin '{filename}': {e}")

    # Load clustering plugins
    if clustering_method != None:
        try:
            folders_clustering = [name for name in os.listdir(plugin_folder + '/clustering') 
                if os.path.isdir(os.path.join(plugin_folder + '/clustering', name))]
            for filename in folders_clustering:
                if clustering_method.lower() == filename.lower():
                    # Load the clustering module
                    clustering_module = importlib.import_module(f'dnabyte.clustering.{filename}.clustering')

                    # Find the class definition in the module
                    for name, obj in inspect.getmembers(clustering_module, inspect.isclass):
                        if issubclass(obj, Cluster) and obj is not Cluster:
                            # Add the class to the plugins dictionary
                            clustering_plugins[filename] = obj
                            break
        except Exception as e:
            print(f"Error loading clustering plugin '{filename}': {e}")

    if recovery_method != None:
        try:
            folders_recovery = [name for name in os.listdir(plugin_folder + '/recovery') 
                if os.path.isdir(os.path.join(plugin_folder + '/recovery', name))]
            for filename in folders_recovery:
                if recovery_method.lower() == filename.lower():
                    # Load the recovery module
                    recovery_module = importlib.import_module(f'dnabyte.recovery.{filename}.recovery')

                    # Find the class definition in the module
                    for name, obj in inspect.getmembers(recovery_module, inspect.isclass):
                        if issubclass(obj, Consensus) and obj is not Consensus:
                            # Add the class to the plugins dictionary
                            recovery_plugins[filename] = obj
                            break
        except Exception as e:
            print(f"Error loading recovery plugin '{filename}': {e}")

    return binarization_plugins, encoding_plugins, storage_plugins, sequencing_plugins, error_plugins, clustering_plugins, recovery_plugins


#TODO: THis needs to be mergerd
def load_synthesis_plugins(synthesis_method):
    """
    Load synthesis plugins from the specified folder and return a dictionary of plugin classes.
    """
    plugin_folder = os.path.dirname(__file__)
    synthesis_plugins = {}  # Define the plugins dictionary
    if isinstance(synthesis_method, str):
        synthesis_method = synthesis_method.lower() if synthesis_method != None else None
    # Normalize method names to lowercase
    if isinstance(synthesis_method, str):
        synthesis_method = synthesis_method.lower() if synthesis_method != None else None

    # Load synthesis plugins
    if synthesis_method != None:
        try:
            folders_synthesize = [name for name in os.listdir(plugin_folder + '/synthesis') 
                if os.path.isdir(os.path.join(plugin_folder + '/synthesis', name))]
            print("folders_synthesize:", folders_synthesize, synthesis_method.lower())
            for filename in folders_synthesize:
                if synthesis_method.lower() == filename.lower():
                    # Load the encode module
                    synthesis_module = importlib.import_module(f'dnabyte.synthesis.{filename}.synthesize')
                    # Find the class definition in the module
                    for name, obj in inspect.getmembers(synthesis_module, inspect.isclass):
                        if issubclass(obj, SimulateSynthesis) and obj is not SimulateSynthesis:
                            # Add the class to the plugins dictionary
                            synthesis_plugins[filename] = obj
                            break
        except Exception as e:
            print(f"Error loading synthesis plugin '{filename}': {e}")

    return synthesis_plugins
