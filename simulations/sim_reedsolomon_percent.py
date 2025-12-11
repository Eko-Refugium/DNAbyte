import os
import random
import string
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.constants import Avogadro

from dnabyte.params import Params
from simulation import Simulation
from auxiliary import create_text_files

# define parameters of the simulation
repeats = 10
reedsoloperc = [y for y in [ 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1] for _ in range(repeats)]


# create files
directory = "./simulations/simfiles"
filenames = create_text_files(directory, [400])

# set other parameters
params = Params.params_range(
        name='synthesis_max_density',
        file_paths=['./simulations/simfiles/textfile_400b.txt'],
        binarization_method='compressed',
        assembly_structure='synthesis',
        encoding_method='max_density',
        mean=100,
        vol=1000000 / Avogadro,
        std_dev=0,
        hybridisation_steps=10000,
        outer_error_correction='reedsolomon',
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        codeword_length=501,
        index_carry_length=34,
        synthesis_method='nosynthpoly',
        sequencing_method='mesa',
        mesa_sequencing_id=68,
        reed_solo_percentage=reedsoloperc,
        seed=69,
        theory='no',
)



# run simulation
sim = Simulation(params)
#res = sim.run(paralel=True, max_workers=8)
res = sim.run(paralel=False)


# extract the success and failure counts for every year
unique_reedsoloperc = list(set(reedsoloperc))
success_failure_counts = {str(reedsoloperc): {'SUCCESS': 0, 'FAILURE': 0} for reedsoloperc in unique_reedsoloperc}

for key, data in res.items():

    index = int(key.split('_')[-1]) - 1

    if res[key].get('status') == 'SUCCESS':
        success_failure_counts[str(params[index].reed_solo_percentage)]['SUCCESS'] += 1
        print('success')
    elif res[key].get('status') == 'FAILURE':
        success_failure_counts[str(params[index].reed_solo_percentage)]['FAILURE'] += 1
        print('failure')

# Print the results
for reedsoloperc, counts in success_failure_counts.items():
    print(f"reedsoloperc: {reedsoloperc}, SUCCESS: {counts['SUCCESS']}, FAILURE: {counts['FAILURE']}")


# Plotting the results
reedsolopercs = list(success_failure_counts.keys())
success_counts = [counts['SUCCESS'] for counts in success_failure_counts.values()]
failure_counts = [counts['FAILURE'] for counts in success_failure_counts.values()]
success_rate = [count / repeats for count in success_counts]

# Sort years and success_rate according to years
sorted_reedsolopercs = sorted([reedsoloperc for reedsoloperc in reedsolopercs])  # Convert years to integers and sort
sorted_success_rate = [success_rate[reedsolopercs.index(str(year))] for year in sorted_reedsolopercs]  # Reorder success_rate


print(success_failure_counts)

print(success_rate)

print(reedsolopercs)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(sorted_reedsolopercs, sorted_success_rate, marker='o', linestyle='-', color='b', label='Success Rate')

# Add labels, title, and legend
# plt.xscale('log')  # Use a logarithmic scale for the x-axis
plt.xlabel('reedsolopercs (log scale)', fontsize=12)
plt.ylabel('Success Rate', fontsize=12)
plt.title('Success Rate vs. reedsolopercs', fontsize=14)
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)


# Set the y-axis range
plt.ylim(0, 1)  # Ensure the y-axis range is always from 0 to 1
plt.xlim(1, max(sorted_reedsolopercs))  # Extend the x-axis limit slightly beyond the maximum year
# Save the plot as an image
job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')  # Generate a unique identifier
output_path = os.path.join('simulations', 'simlogs', f'sim_storage_{job_identifier}.png')
plt.tight_layout()
plt.savefig(output_path, dpi=300)  # Save the plot with high resolution
print(f"Plot saved as {output_path}")

# Show the plot
plt.show()