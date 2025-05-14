import os
import random
import string
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.constants import Avogadro
import pickle

from dnabyte.params import Params
from simulations.simulation import Simulation
from simulations.auxiliary import create_text_files


# define parameters of the simulation
repeats = 10
iid_error_rate = [y for y in [0] for _ in range(repeats)]


# set other parameters
params = Params.params_range(
        name='sim_linear_chain_noEC',
        filename='./simulations/simfiles/textfile_40b_restored.txt',
        assembly_structure='linear_assembly',
        encoding_scheme='linear_encoding',
        library_name='lib_simple_30_400.csv',
        mean=10,
        vol=1000000 / Avogadro,
        std_dev=0,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction='reedsolomon',
        dna_barcode_length=5,  
        codeword_maxlength_positions=5,
        years=0,
        storage_conditions='biogene',
        codeword_length=100,
        percent_of_symbols=2,
        ltcode_header=4,
        index_carry_length=3,
        synthesis_method=None,
        sequencing_method='iid',
        iid_error_rate=iid_error_rate,
        reed_solo_percentage=0.9,
        sigma_amount=None,
        seed=42,
        theory='no'
)
filename='./simulations/simfiles/Order_lib_simple_30_144_1.txt'
encoding_scheme='linear_encoding'
mean=5

filename='./simulations/simfiles/textfile_40b.txt'
encoding_scheme='linear_encoding'
mean=5

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
    file.write("iid_error_rate: " + str(iid_error_rate) +", success_rate: " + str(success_rate) + ", success_failure_counts: " + str(success_failure_counts) + ", library" + str(filename) + ", encodingscheme" + str(encoding_scheme) + ", mean" + str(mean) + "\n")

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