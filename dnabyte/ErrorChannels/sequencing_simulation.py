import os
import json
import random
import numpy as np
from flask import current_app

from dnabyte.ErrorChannels.simulator_api import Graph
from dnabyte.ErrorChannels.sequencing_error import SequencingError

def remove_spaces_from_nested_list(nested_list):
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


def process_sequences(sequences, err_att_syn, err_rate_syn):
    
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
            seed, error_dict = sequencing_err.lit_error_rate_mutations()
            mod_seq = g.graph.nodes[0]['seq']
            mod_seq = mod_seq.replace(' ', '')

            return mod_seq, error_dict

    modified_sequences, error_dicts = process_element(sequences)
    return modified_sequences, error_dicts



def sequencing_simulation(sequences, method_id, errorrate=None):
    
    if method_id == 'iid':
        sequenceserror = sequences
        for i in range(len(sequenceserror)):
            for j in range(len(sequenceserror[i])):
                for k in range(len(sequenceserror[i][j])):
                    if random.random() < errorrate:
                        # randomly select a new base
                        new_base = random.choice(['A', 'C', 'G', 'T'])
                        # replace the base at the mutation position with the new base
                        sequenceserror[i][j] = (sequences[i][j][:k] + new_base + sequences[i][j][k + 1:])      
               
                
        error_dict = {}
        return sequenceserror,error_dict
    
    else:
        #flask version
        # parent_dir = os.path.dirname(current_app.root_path)
        # nonflask
        parent_dir = os.getcwd()
        
        # Get the absolute path of the file relative to the parent directory
        file_path = os.path.join(parent_dir, 'mi_dna_disc', 'ErrorChannels', 'seq_table.json')

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
        sequences_modified, error_dict =  process_sequences(sequences, err_att_syn, err_rate_syn)
            
        return sequences_modified, error_dict
