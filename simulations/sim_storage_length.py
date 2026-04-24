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
repeats = 100
years = [y for y in [1, 10, 100, 1000, 10000] for _ in range(repeats)]

# create files
directory = "./simulations/simfiles"
filenames = create_text_files(directory, [400])

# set other parameters
params = Params.params_range(
        name='synthesis_max_density',
        file_paths=['./simulations/simfiles/textfile_400b.txt'],
        assembly_structure='synthesis',
        encoding_method='max_density',
        binarization_method='compressed',
        mean=100,
        vol=1000000 / Avogadro,
        std_dev=0,
        hybridisation_steps=10000,
        inner_error_correction=None,
        outer_error_correction='reedsolomon',
        dna_barcode_length=34,  
        codeword_maxlength_positions=18,
        years=years,
        storage_conditions='biogene',
        codeword_length=501,
        percent_of_symbols=6,
        index_carry_length=34,
        # synthesis_method='mesa',
        # mesa_synthesis_id=68,
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


# extract the success and failure counts for every year
unique_years = list(set(years))
success_failure_counts = {str(year): {'SUCCESS': 0, 'FAILURE': 0} for year in unique_years}

for key, data in res.items():
    index = int(key.split('_')[-1]) - 1
    if res[key].get('status') == 'SUCCESS':
        success_failure_counts[str(params[index].years)]['SUCCESS'] += 1
    elif res[key].get('status') == 'FAILURE':
        success_failure_counts[str(params[index].years)]['FAILURE'] += 1

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

# Set the y-axis range
plt.ylim(0, 1)  # Ensure the y-axis range is always from 0 to 1
plt.xlim(1, max(sorted_years))  # Extend the x-axis limit slightly beyond the maximum year
# Save the plot as an image
job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')  # Generate a unique identifier
output_path = os.path.join('simulations', 'simlogs', f'sim_storage_{job_identifier}.png')
plt.tight_layout()
plt.savefig(output_path, dpi=300)  # Save the plot with high resolution
print(f"Plot saved as {output_path}")

# Show the plot
plt.show()