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
years = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192] 
directory = "./simulations/simfiles"
filenames = create_text_files(directory, [40])

# set remaining parameters
params = Params.params_range(
        name='success_rate_of_synthesis_max_density',
        filename='./simulations/simfiles/textfile_40b.txt',
        assembly_structure='synthesis',
        encoding_scheme='max_density_encoding',
        library_name='',
        mean=200,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction='ltcode',
        outer_error_correction='reedsolomon',
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        years=years,
        storage_conditions=None,
        codeword_length=501,
        percent_of_symbols=2,
        index_carry_length=34,
        synthesis_method=68,
        ltcode_header=34,
        sequencing_method=None,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        theory='no'
)

print(params)

# run simulation
sim = Simulation(params, debug=True)
res = sim.run()

# delete files
for filename in filenames:
    os.remove(filename)
    #os.remove(os.path.join(directory, filename))

# Sum up the durations of every step for each element of res
success_rate = []

for name, steps in res.items():
    total_duration = sum(step_info['duration'] for step_info in steps.values() if isinstance(step_info, dict) and 'duration' in step_info)
    success_rate.append(total_duration)

print(res)

# # Create the plot
# plt.figure(figsize=(10, 6))
# plt.plot(years, total_durations, marker='o')
# plt.xlabel('File Size (bytes)')
# plt.ylabel('Total Duration (seconds)')
# plt.title('File Size vs. Total Duration')
# plt.grid(True)

# # Save the plot to a file
# plt.savefig(os.path.join(directory, 'sim_duration_plot.png'))

# # Show the plot
# plt.show()