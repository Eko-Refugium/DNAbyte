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
years = [y for y in [10, 100, 1000, 10000, 20000, 40000, 60000, 80000, 90000, 100000, 110000, 120000] for _ in range(repeats)]


# set other parameters
params = Params.params_range(
        name='synthesis_max_density',
        file_paths=['./simulations/simfiles/textfile_40b.txt'],
        assembly_structure='synthesis',
        encoding_method='max_density',
        mean=200,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction='reedsolomon',
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        binarization_method='compressed',
        years=years,
        storage_conditions='biogene',
        codeword_length=501,
        percent_of_symbols=2,
        index_carry_length=34,
        # synthesis_method='mesa',
        # mesa_synthesis_id=68,
        # sequencing_method='mesa',
        # mesa_sequencing_id=41,
        reed_solo_percentage=0.9,
        seed=42,
        theory='no'
)

# run simulation
sim = Simulation(params)
res = sim.run(paralel=False)

# Debug: Check the structure of res
if len(res) > 0:
    first_key = list(res.keys())[0]
    print(f"Sample key: {first_key}")
    print(f"Sample value: {res[first_key]}")
    print(f"Total results: {len(res)}")

# extract the success and failure counts for every year
unique_years = list(set(years))
success_failure_counts = {str(year): {'SUCCESS': 0, 'FAILURE': 0} for year in unique_years}

for key, data in res.items():

    index = int(key.split('_')[-1]) - 1

    if res[key].get('status') == 'SUCCESS':
        success_failure_counts[str(params[index].years)]['SUCCESS'] += 1
        print('success')
    elif res[key].get('status') == 'FAILURE':
        success_failure_counts[str(params[index].years)]['FAILURE'] += 1
        print('failure')

# Print the results
for year, counts in success_failure_counts.items():
    print(f"Year: {year}, SUCCESS: {counts['SUCCESS']}, FAILURE: {counts['FAILURE']}")


# Plotting the results
years = list(success_failure_counts.keys())
success_counts = [counts['SUCCESS'] for counts in success_failure_counts.values()]
failure_counts = [counts['FAILURE'] for counts in success_failure_counts.values()]
success_rate = [count / repeats for count in success_counts]

# Sort years and success_rate according to years
sorted_years = sorted([int(year) for year in years])  # Convert years to integers and sort
sorted_success_rate = [success_rate[years.index(str(year))] for year in sorted_years]  # Reorder success_rate


print(success_failure_counts)

print(success_rate)

print(years)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(sorted_years, sorted_success_rate, marker='o', linestyle='-', color='b', label='Success Rate')

# Add labels, title, and legend
plt.xscale('log')  # Use a logarithmic scale for the x-axis
plt.xlabel('Years (log scale)', fontsize=12)
plt.ylabel('Success Rate', fontsize=12)
plt.title('Success Rate vs. Years', fontsize=14)
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