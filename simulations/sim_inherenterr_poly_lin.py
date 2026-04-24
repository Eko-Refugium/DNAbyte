import os
import random
import string
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.constants import Avogadro
import pickle

from dnabyte.params import Params
from simulation import Simulation
from auxiliary import create_text_files


# define parameters of the simulation
repeats = 10
iid_error_rate = [y for y in [0,0.02,0.04,0.06,0.08,0.1] for _ in range(repeats)]


# set other parameters
params = Params.params_range(
        name='sim_linear_chain_noEC',
        file_paths=['./simulations/simfiles/textfile_40b.txt'],
        binarization_method='compressed',
        assembly_structure='positional_assembly',
        encoding_method='poly_chain',
        library_name='polymeraselibinparts_500(40)_100(35)_70.csv',
        mean=1,
        std_dev=0,
        hybridisation_steps=10000,
        dna_barcode_length=2,  
        codeword_maxlength_positions=2,
        codeword_length=100,
        sequencing_method='iid',
        iid_error_rate=iid_error_rate,
        seed=42,
        theory='no'
)

filename='./simulations/simfiles/textfile_640b.txt'
sim = Simulation(params)
#res = sim.run(paralel=True, max_workers=8)
res = sim.run(paralel=False)

# extract the success and failure counts for every year
unique_iid_error_rates = list(set(iid_error_rate))
success_failure_counts = {str(iid_error_rate): {'SUCCESS': 0, 'FAILURE': 0} for iid_error_rate in unique_iid_error_rates}

for key, data in res.items():

    index = int(key.split('_')[-1]) - 1

    if res[key].get('status') == 'SUCCESS':
        success_failure_counts[str(params[index].iid_error_rate)]['SUCCESS'] += 1
        print('success')
    elif res[key].get('status') == 'FAILURE':
        success_failure_counts[str(params[index].iid_error_rate)]['FAILURE'] += 1
        print('failure')

# Print the results
for iid_error_rate, counts in success_failure_counts.items():
    print(f"iid_error_rate: {iid_error_rate}, SUCCESS: {counts['SUCCESS']}, FAILURE: {counts['FAILURE']}")


# Plotting the results
iid_error_rate = list(success_failure_counts.keys())
success_counts = [counts['SUCCESS'] for counts in success_failure_counts.values()]
failure_counts = [counts['FAILURE'] for counts in success_failure_counts.values()]
success_rate = [count / repeats for count in success_counts]

# Sort years and success_rate according to years
sorted_iid_error_rates = sorted([iid_error_rate for iid_error_rate in iid_error_rate])  # Convert years to integers and sort
sorted_success_rate = [success_rate[iid_error_rate.index(str(year))] for year in sorted_iid_error_rates]  # Reorder success_rate


print(success_failure_counts)

print(success_rate)

print(iid_error_rate)

with open("dataofsims.txt", "a+") as file:
    file.write("iid_error_rate: " + str(iid_error_rate) +", success_rate: " + str(success_rate) + ", success_failure_counts: " + str(success_failure_counts) + ", library" + filename + "\n")

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(sorted_iid_error_rates, sorted_success_rate, marker='o', linestyle='-', color='b', label='Success Rate')

# Add labels, title, and legend
# plt.xscale('log')  # Use a logarithmic scale for the x-axis
plt.xlabel('iid_error_rates', fontsize=12)
plt.ylabel('Success Rate', fontsize=12)
plt.title('Success Rate vs. iid_error_rates', fontsize=14)
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Save the plot as an image
job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')  # Generate a unique identifier
output_path = os.path.join('simulations', 'simlogs', f'sim_storage_{job_identifier}.png')
plt.tight_layout()
plt.savefig(output_path, dpi=300)  # Save the plot with high resolution
print(f"Plot saved as {output_path}")

# Show the plot
plt.show()