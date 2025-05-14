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
repeats = 100
years = [y for y in [10, 100, 1000, 10000, 20000, 40000, 60000, 80000, 90000, 100000, 110000, 120000] for _ in range(repeats)]


# set other parameters
params = Params.params_range(
        name='synthesis_max_density',
        filename='./simulations/simfiles/textfile_40b.txt',
        assembly_structure='synthesis',
        encoding_scheme='max_density_encoding',
        library_name='',
        mean=200,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction='ltcode',
        outer_error_correction='reed_solomon',
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        years=years,
        storage_conditions='biogene',
        codeword_length=501,
        percent_of_symbols=2,
        index_carry_length=34,
        synthesis_method=68,
        sequencing_method=None,
        reed_solo_percentage=0.9,
        sigma_amount=None,
        seed=42,
        theory='no'
)

job_identifier = '20250407_084254'

pickle_file_path = os.path.join('simulations', 'simlogs', f'res_{job_identifier}.pickle')

# Load the res object
try:
    with open(pickle_file_path, 'rb') as pickle_file:
        res = pickle.load(pickle_file)
    print(f"Loaded res object from {pickle_file_path}")
except FileNotFoundError:
    print(f"Pickle file not found: {pickle_file_path}")
except Exception as e:
    print(f"An error occurred while loading the pickle file: {e}")

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