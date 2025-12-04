import os
import json
import numpy as np

from dnabyte.synthesize import SimulateSynthesis
from .error_graph import Graph
from .sequencing_error import SequencingError


class MESA(SimulateSynthesis):
    """
    MESA (Molecular Error Simulation Algorithm) is a class that simulates sequencing errors in DNA sequences.
    It uses a graph-based approach to model the sequencing process and introduces errors based on specified parameters.
    """

    def simulate(self, data):
        """
        Simulate sequencing errors using the MESA model.
        
        :param data: A list of DNA sequences.
        :return: A list of sequenced DNA sequences.
        """
        parent_dir = os.getcwd()
        
        method_id = self.params.mesa_synthesis_id

        # Get the absolute path of the file relative to the parent directory
        file_path = os.path.join(parent_dir, 'dnabyte', 'synthesis', 'mesa',  'syn_table.json')

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
        sequences_multiples = [seq for seq in data for _ in range(max(1, int(np.random.normal(self.params.mean, self.params.std_dev))))]

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

        average_copy_number = len(sequences_modified) / len(data)
        info = {
            "average_copy_number": average_copy_number,
            "number_of_synthesis_errors": "unknown",
            "error_dict": {}
        }

        return sequences_modified, info
        
def attributes(params):

    if 'mean' not in params.__dict__ or params.mean is None:
        mean = 10
    else:
        mean = params.mean

    if 'std_dev' not in params.__dict__ or params.std_dev is None:
        std_dev = 0
    else:
        std_dev = params.std_dev

    if 'mesa_synthesis_id' not in params.__dict__ or params.mesa_synthesis_id is None:
        mesa_synthesis_id = 68
    else:
        if params.mesa_synthesis_id not in [3, 4, 5, 6, 7, 68, 69, 70, 71, None]:
            raise ValueError("Invalid Synthesis ID Error")
        else:
            mesa_synthesis_id = params.mesa_synthesis_id

    return {"mean": mean, 
            "std_dev": std_dev, 
            "mesa_synthesis_id": mesa_synthesis_id}