import os
import importlib
import inspect

from dnabyte.encode import Encode
from dnabyte.store import SimulateStorage

from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.synthesize import SimulateSynthesis


def load_plugins(encoding_method, synthesis_method, storage_conditions, sequencing_method, binarization_method):
        """
        Load plugins from the specified folder and return a dictionary of plugin classes.
        """
        plugin_folder = os.path.dirname(__file__)
        encoding_plugins = {}  # Define the plugins dictionary
        synthesis_plugins = {}  # Define the plugins dictionary
        storage_plugins = {}  # Define the plugins dictionary
        sequencing_plugins = {}  # Define the plugins dictionary
        binarization_plugins = {}  # Define the plugins dictionary

        encoding_method = encoding_method.lower() if encoding_method != None else None
        if isinstance(synthesis_method, str):
            synthesis_method = synthesis_method.lower() if synthesis_method != None else None

        if isinstance(storage_conditions, str):
            storage_conditions = storage_conditions.lower() if storage_conditions != None else None
        elif isinstance(storage_conditions, list):
            storage_conditions = [condition.lower() for condition in storage_conditions]
        sequencing_method = sequencing_method.lower() if sequencing_method != None else None
        binarization_method = binarization_method.lower() if binarization_method != None else None
        path_encode = plugin_folder + '/encoding'

        folders_encode = [name for name in os.listdir(path_encode) 
            if os.path.isdir(os.path.join(path_encode, name))]
        

        # Load encoding plugins
        if encoding_method != None:
            try:
                for filename in folders_encode:
                    # if filename.startswith('encode_') and filename.endswith('.py'):
                    # encoding_type = filename[7:-3]  # Extract the encoding type (e.g., 'maxdensity' from 'encode_maxdensity.py')
                    print(filename.lower(), encoding_method, 'filename in load plugins')
                    # Load the encode module
                    if encoding_method.lower() == filename.lower():
                        #print('dnabyte.EncodingMethods' + f'.{filename}.encode')
                        encode_module = importlib.import_module('dnabyte.encoding' + f'.{filename}.encode')
                        #print(encode_module)
                        # Find the class definition in the module
                        # encode_module.validate()
                        for name, obj in inspect.getmembers(encode_module, inspect.isclass):
                            if issubclass(obj, Encode) and obj is not Encode:
                                # Add the class to the plugins dictionary
                                encoding_plugins[filename] = obj
                                break
            except Exception as e:
                print(f"Error loading encoding plugin '{filename}': {e}")
        # else:
        #     raise ValueError("No encoding method specified.")

        if sequencing_method != None:
            try:
                for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/sequencing_plugins') if f.endswith('.py')]:
                    # if filename.startswith('encode_') and filename.endswith('.py'):
                    module_name = filename[:-3] 
                    if sequencing_method.lower() == module_name.lower():
                        # Load the encode module
                        sequencing_module = importlib.import_module('dnabyte.sequencing_plugins.synthesis_' + f'{module_name}')
                        # Find the class definition in the module
                        for name, obj in inspect.getmembers(sequencing_module, inspect.isclass):
                            if issubclass(obj, SimulateSequencing) and obj is not SimulateSequencing:
                                # Add the class to the plugins dictionary
                                sequencing_plugins[module_name] = obj
                                break
            except Exception as e:
                print(f"Error loading sequencing plugin '{module_name}': {e}")
        else:
            print("No sequencing method specified, skipping sequencing plugin loading.")

        if synthesis_method != None:
            try:
                for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/synthesis') if f.endswith('.py')]:
                    # if filename.startswith('encode_') and filename.endswith('.py'):
                    module_name = filename[:-3] 
                    if synthesis_method.lower() == module_name.lower():
                        # Load the encode module
                        synthesis_module = importlib.import_module('dnabyte.synthesis.synthesis_' + f'{module_name}')
                        # Find the class definition in the module
                        for name, obj in inspect.getmembers(synthesis_module, inspect.isclass):
                            if issubclass(obj, SimulateSynthesis) and obj is not SimulateSynthesis:
                                # Add the class to the plugins dictionary
                                synthesis_plugins[module_name] = obj
                                break
            except Exception as e:
                print(f"Error loading synthesis plugin '{module_name}': {e}")
        else:
            print("No synthesis method specified, skipping synthesis plugin loading.")

        print('STORAGE CONDITIONS:', storage_conditions)

        if storage_conditions != None:
            try:   
                for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/storage') if f.endswith('.py')]:
                    # if filename.startswith('encode_') and filename.endswith('.py'):
                    print('filename:', filename)
                    module_name = filename[:-3]
                    print('module_name:', module_name)

                    if isinstance(storage_conditions, str):
                        if storage_conditions.lower() == module_name.lower():
                            # Load the encode module
                            storage_module = importlib.import_module('dnabyte.storage.storage_' + f'{module_name}')
                            # Find the class definition in the module
                            for name, obj in inspect.getmembers(storage_module, inspect.isclass):
                                if issubclass(obj, SimulateStorage) and obj is not SimulateStorage:
                                    # Add the class to the plugins dictionary
                                    storage_plugins[module_name] = obj
                                    break
                    elif isinstance(storage_conditions, list):
                        for condition in storage_conditions:
                            if 'storage_' + condition.lower() == module_name.lower():
                                # Load the encode module
                                storage_module = importlib.import_module('dnabyte.storage.' + f'{module_name}')
                                # Find the class definition in the module
                                for name, obj in inspect.getmembers(storage_module, inspect.isclass):
                                    if issubclass(obj, SimulateStorage) and obj is not SimulateStorage:
                                        # Add the class to the plugins dictionary
                                        storage_plugins[module_name] = obj
                                        break
            except Exception as e:
                print(f"Error loading storage plugin '{module_name}': {e}")
        else:
            print("No storage conditions specified, skipping storage plugin loading.")

        if binarization_method != None:
            try:
                for filename in [f for f in os.listdir(os.path.dirname(__file__) + '/binarization') if f.endswith('.py')]:
                    #print(filename)
                    module_name = filename[:-3] 
                    if binarization_method.lower() == module_name.lower():
                        
                        # Load the encode module
                        binarization_module = importlib.import_module('dnabyte.binarization.' + f'{module_name}')
                        # Find the class definition in the module
                        for name, obj in inspect.getmembers(binarization_module, inspect.isclass):
                            if issubclass(obj, Binarize) and obj is not Binarize:
                            # Add the class to the plugins dictionary
                                binarization_plugins[module_name] = obj
                                break
            except Exception as e:
                print(f"Error loading binarization plugin '{module_name}': {e}")
        # else:
        #     raise ValueError("No binarization method specified.")


        return encoding_plugins, synthesis_plugins, storage_plugins, sequencing_plugins, binarization_plugins

