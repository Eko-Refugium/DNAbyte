import os
import pickle
import matplotlib.pyplot as plt

# Define the pickle file names and corresponding labels
pickle_files = [
    ('res_20250407_170104.pickle', 'res_6'),
    ('res_20250407_165941.pickle', 'res_5'),
    ('res_20250407_165833.pickle', 'res_4'),
    ('res_20250407_165728.pickle', 'res_3'),
    ('res_20250407_165626.pickle', 'res_2'),
    ('res_20250407_165437.pickle', 'res_1')
]

# Initialize a dictionary to store results
results = {}

# Load each pickle file
for file_name, label in pickle_files:
    try:
        with open(os.path.join('simulations', 'simlogs', file_name), 'rb') as pickle_file:
            results[label] = pickle.load(pickle_file)
        print(f"Loaded {label} from {file_name}")
    except Exception as e:
        print(f"Failed to load {file_name}: {e}")

# Define parameters
repeats = 200
years = [y for y in [1, 10, 100, 1000, 10000, 20000, 40000, 60000, 80000, 90000, 100000, 110000, 120000] for _ in range(repeats)]

# Plot all curves
plt.figure(figsize=(10, 6))

for label, res in results.items():
    # Extract success and failure counts
    unique_years = list(set(years))
    success_failure_counts = {str(year): {'SUCCESS': 0, 'FAILURE': 0} for year in unique_years}

    for key, data in res.items():
        index = int(key.split('_')[-1]) - 1

        if res[key].get('status') == 'SUCCESS':
            success_failure_counts[str(years[index])]['SUCCESS'] += 1
        elif res[key].get('status') == 'FAILURE':
            success_failure_counts[str(years[index])]['FAILURE'] += 1

    # Calculate success rates
    years_list = list(success_failure_counts.keys())
    success_counts = [counts['SUCCESS'] for counts in success_failure_counts.values()]
    success_rate = [count / repeats for count in success_counts]

    # Sort years and success_rate according to years
    sorted_years = sorted([int(year) for year in years_list])  # Convert years to integers and sort
    sorted_success_rate = [success_rate[years_list.index(str(year))] for year in sorted_years]  # Reorder success_rate

    # Plot the curve
    plt.plot(sorted_years, sorted_success_rate, marker='o', linestyle='-', label=label)

# Add labels, title, and legend
plt.xscale('log')  # Use a logarithmic scale for the x-axis
plt.xlabel('Years (log scale)', fontsize=12)
plt.ylabel('Success Rate', fontsize=12)
plt.title('Success Rates for Increasing values of the LT Code Parameter', fontsize=14)
plt.legend()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

# Set the y-axis range
plt.ylim(0, 1)  # Ensure the y-axis range is always from 0 to 1
plt.xlim(1, max(sorted_years))  # Extend the x-axis limit slightly beyond the maximum year

# Save the plot as an image
output_path = os.path.join('simulations', 'simlogs', 'combined_sim_storage.png')
plt.tight_layout()
plt.savefig(output_path, dpi=300)  # Save the plot with high resolution
print(f"Combined plot saved as {output_path}")

# Show the plot
plt.show()