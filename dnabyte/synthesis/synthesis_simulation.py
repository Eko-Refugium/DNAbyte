
import os
import json
import random
import numpy as np

from flask import current_app
from dnabyte.sequencing.simulator_api import Graph
from dnabyte.sequencing.sequencing_error import SequencingError

def synthesis_simulation(sequences, method_id, mean, std_dev):

    # Get the parent directory of the current app root path,
    #flask version
    # parent_dir = os.path.dirname(current_app.root_path)
    
    # nonflask
    parent_dir = os.getcwd()
    
    # Get the absolute path of the file relative to the parent directory
    file_path = os.path.join(parent_dir, 'dnabyte', 'ErrorChannels', 'syn_table.json')

    # get the error parameters for the designated method
    # get dictionary of synthesis parameters from JSON file
    with open(file_path, 'r') as f:
        synth_dict = json.load(f)

    # Subset the dictionary
    method_dict = {item['id']: item for item in synth_dict if item.get('id') == method_id}       

    # get the error parameters for the designated method
    err_rate_syn = method_dict[method_id]['err_data']
    err_att_syn = method_dict[method_id]['err_attributes']
    method_id = method_dict[method_id].get('id')

    # create copies of sequences based on normal distribution
    sequences_multiples = [seq for seq in sequences for _ in range(max(1, int(np.random.normal(mean, std_dev))))]

    sequences_modified = []

    for original_sequence in sequences_multiples:

        # generate seed
        org_seed = int(np.random.randint(0, 4294967295, dtype=np.uint32))
        seed = np.uint32(float(org_seed) % 4294967296) if org_seed else None

        # The Graph for all types of errors
        g = Graph(None, original_sequence)

        synth_err = SequencingError(original_sequence, g, 'synthesis', err_att_syn, err_rate_syn, seed=seed)
        seed = synth_err.lit_error_rate_mutations()

        mod_seq = g.graph.nodes[0]['seq']
        mod_seq = mod_seq.replace(' ', '')

        sequences_modified.append(mod_seq)

    return sequences_modified

