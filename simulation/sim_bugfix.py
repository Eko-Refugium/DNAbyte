import os
import random
import string
import matplotlib.pyplot as plt
from scipy.constants import Avogadro

from dnabyte.params import Params
from simulations.simulation import Simulation
from simulations.auxiliary import create_text_files

# define parameters of the simulation

# create files
sizes = [2560]  # Sizes: 10, 20, 40, 80, ..., 40960 bytes
#sizes = [40]  # Sizes: 10, 20, 40, 80, ..., 40960 bytes

directory = "./simulations/simfiles"
filenames = create_text_files(directory, sizes)

# set remaining parameters
params = Params.params_range(
        name='synthesis_max_density',
        filename=filenames,
        assembly_structure='synthesis',
        encoding_scheme='max_density_encoding',
        library_name='',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction='ltcode',
        outer_error_correction='reed_solomon',
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        years=0,
        storage_conditions=None,
        codeword_length=501,
        percent_of_symbols=2,
        index_carry_length=34,
        synthesis_method=68,
        sequencing_method=None,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        theory='no'
)

# run simulation
sim = Simulation(params)
res = sim.run()

# delete files
for filename in filenames:
    os.remove(filename)
    #os.remove(os.path.join(directory, filename))

