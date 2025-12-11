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
sizes = [40 * (2 ** i) for i in range(17)]  # Sizes: 10, 20, 40, 80, ..., 40960 bytes
directory = "./simulations/simfiles"
filenames = create_text_files(directory, sizes)


# set remaining parameters
params = Params.params_range(
        name='Linear_chain',
        file_paths=filenames,
        assembly_structure='linear_assembly',
        encoding_method='linear_encoding',
        binarization_method='compressed',
        library_name='20bp_Lib.csv',
        inner_error_correction=None,
        outer_error_correction='reedsolomon',
        reed_solo_percentage=0.8,
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        codeword_length=501

)

# run simulation
sim = Simulation(params)
res = sim.run()

# delete files
for filename in filenames:
    os.remove(filename)
    #os.remove(os.path.join(directory, filename))

# Sum up the durations of every step for each element of res
total_durations = []

for name, steps in res.items():
    total_duration = sum(step_info['duration'] for step_info in steps.values() if isinstance(step_info, dict) and 'duration' in step_info)
    total_durations.append(total_duration)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(sizes, total_durations, marker='o')
plt.xlabel('File Size (bytes)')
plt.ylabel('Total Duration (seconds)')
plt.title('File Size vs. Total Duration')
plt.grid(True)

# Save the plot to a file
plt.savefig(os.path.join(directory, 'sim_duration_ec_plot.png'))

# Show the plot
plt.show()