"""Simulation for de Bruijn graph recovery depth."""

#Function that makes a random n bit txt file
import os
import string
import random

from simulations.simulation import Simulation


def create_random_file(output_path: str, size_bytes: int = 500):
    """Create a random text file for encoding"""
    content = ''.join(random.choices(string.ascii_letters + string.digits + ' \n', k=size_bytes))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(content)
    return output_path

# run one simulation to check if it works
def one_sim(encoding, params):
    # Run the encoding and decoding simulation for the given encoding method
    # This is a placeholder for the actual simulation logic
    # Replace this with the actual implementation of your encoding and decoding process
    try:
        sim = Simulation(params)
        run = sim.run(paralel=False)
        for key, data in run.items():
            if run[key].get('status') == 'SUCCESS':
                print(f"Simulation for {encoding} succeeded.")
                return True
            elif run[key].get('status') == 'FAILURE':
                raise RuntimeError(f"Simulation for {encoding} failed with error: {run[key].get('error')}")
    except Exception as e:
        print(f"Error during simulation for {encoding}: {e}")
        return False



# run a big simulation with errors of diffrent magnitudes 


# 