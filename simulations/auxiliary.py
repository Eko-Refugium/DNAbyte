import random
import string
import os

def generate_random_ascii(size):
    """
    Generates a random ASCII string of the given size.
    
    Parameters:
    size (int): The size of the ASCII string to generate.
    
    Returns:
    str: A random ASCII string of the given size.
    """
    return ''.join(random.choices(string.ascii_letters + string.digits + string.punctuation, k=size))

def create_text_files(directory, sizes):
    """
    Creates text files with random ASCII characters of sizes 10 bytes, 20 bytes, 40 bytes, ..., 40960 bytes.
    
    Parameters:
    directory (str): The directory where the text files will be created.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    filenames = []

    for size in sizes:
        filename = os.path.join(directory, f"textfile_{size}b.txt")
        filenames.append(filename)
        with open(filename, 'w') as f:
            f.write(generate_random_ascii(size))

    return filenames

def create_one_text_file(directory, size):
    """
    Creates text files with random ASCII characters of sizes 10 bytes, 20 bytes, 40 bytes, ..., 40960 bytes.
    
    Parameters:
    directory (str): The directory where the text files will be created.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    

    filename = f"textfile_{size}b.txt"
    print(directory+filename)
    with open(directory+filename, 'w') as f:
        f.write(generate_random_ascii(size))
    print(f"Created {filename} with size {size} bytes")

    return filename