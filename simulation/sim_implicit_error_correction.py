import os
import random
import string
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.constants import Avogadro

from dnabyte.params import Params
from simulations.simulation import Simulation
from simulations.auxiliary import create_text_files

# define parameters of the simulation
repeats = 100
years = [1]


# create files
directory = "./simulations/simfiles"
filenames = create_text_files(directory, [400])


# set other parameters
params = Params.params_range(
        name='synthesis_max_density',
        filename='./simulations/simfiles/textfile_400b.txt',
        assembly_structure='linear_assembly',
        encoding_scheme='linear_encoding',
        library_name='20bp_Lib.csv',
        mean=200,
        vol=1000000 / Avogadro,
        std_dev=0,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction=None,
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        years=1,
        storage_conditions=None,
        codeword_length=501,
        percent_of_symbols=6,
        index_carry_length=34,
        synthesis_method=None,
        sequencing_method=None,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        seed=69,
        theory='no',
        ltcode_header=34
)

# run simulation
sim = Simulation(params)
res = sim.run(paralel=False)

print(res)