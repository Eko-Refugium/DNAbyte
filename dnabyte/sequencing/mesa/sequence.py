import os
import json
import numpy as np

from dnabyte.sequence import SimulateSequencing
from .simulator import Graph
from .sequencing_error import SequencingError

class MESA(SimulateSequencing):
    """
    Simulate sequencing errors using the MESA model.
    This model is based on the MESA (Molecular Error Simulation Algorithm) approach.
    """

    def simulate(self, data):
        """
        Simulate sequencing errors using the MESA model.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        # TODO: the sequencing simulator should not be doing coverage simulation
        # Create copies of sequences based on normal distribution (coverage)
        #sequences_multiples = [seq for seq in data for _ in range(max(1, int(np.random.normal(self.params.mean, self.params.std_dev))))]
        
        # Call the sequencing simulation function
        modified_sequences, error_dicts = self.sequencing_simulation(data, method_id=self.params.mesa_sequencing_id)
                
        info = {
            "average_copy_number": 1.0,  # No coverage simulation in sequencing
            "number_of_sequencing_errors": sum([len(d) for d in error_dicts]),
            "error_dict": error_dicts
        }
        
        return modified_sequences, info

    #TODO: is this function needed?
    def remove_spaces_from_nested_list(self, nested_list):
        def remove_spaces(s):
            return s.replace(' ', '')

        def process_list(lst):
            for i in range(len(lst)):
                if isinstance(lst[i], list):
                    process_list(lst[i])
                elif isinstance(lst[i], str):
                    lst[i] = remove_spaces(lst[i])
        process_list(nested_list)
        return nested_list


    def process_sequences(self, sequences, err_att_syn, err_rate_syn):
        
        def process_element(element):
            if isinstance(element, list):
                results = [process_element(sub_element) for sub_element in element]
                mod_seqs, error_dicts = zip(*results)
                return list(mod_seqs), list(error_dicts)
            elif isinstance(element, str):
                # generate seed
                org_seed = int(np.random.randint(0, 4294967295, dtype=np.uint32))
                seed = np.uint32(float(org_seed) % 4294967296) if org_seed else None

                # The Graph for all types of errors
                g = Graph(None, element)

                sequencing_err = SequencingError(element, g, 'sequencing', err_att_syn, err_rate_syn, seed=seed)
                seed = sequencing_err.lit_error_rate_mutations()
                mod_seq = g.graph.nodes[0]['seq']
                mod_seq = mod_seq.replace(' ', '')

                # Create error dict (currently empty, could be populated from graph if needed)
                error_dict = {}

                return mod_seq, error_dict

        modified_sequences, error_dicts = process_element(sequences)
        return modified_sequences, error_dicts


    def sequencing_simulation(self, sequences, method_id):

        # Get the directory where this module is located
        module_dir = os.path.dirname(os.path.abspath(__file__))

        # Get the absolute path of the JSON file in the same directory as this module
        file_path = os.path.join(module_dir, 'seq_table.json')

        # get the error parameters for the designated method
        # get dictionary of synthesis parameters from JSON file
        
        with open(file_path, 'r') as f:
            seq_dict = json.load(f)

        # Subset the dictionary
        method_dict = {item['id']: item for item in seq_dict if item.get('id') == method_id}       

        # get the error parameters for the designated method
        err_rate_syn = method_dict[method_id]['err_data']
        err_att_syn = method_dict[method_id]['err_attributes']
        method_id = method_dict[method_id].get('id')
        sequences_modified, error_dict =  self.process_sequences(sequences, err_att_syn, err_rate_syn)
            
        return sequences_modified, error_dict

def attributes(params):
    # required
    if 'mesa_sequencing_id' not in params.__dict__ or params.mesa_sequencing_id is None:
        # set default sequencing method id
        mesa_sequencing_id = 38
    else:
        if params.mesa_sequencing_id not in [41, 40, 37, 36, 39, 38, 35]:
            raise ValueError("Invalid sequencing method id")
        else:
            mesa_sequencing_id = params.mesa_sequencing_id
        
    return {"mesa_sequencing_id": mesa_sequencing_id}