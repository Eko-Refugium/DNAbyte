from scipy.constants import Avogadro
from dnabyte.params import Params
from dnabyte.library import Library

simulation_parameters = [
    Params(
        name='synthesis_max_density',
        filename=[
            'textfile_10b.txt',
            'textfile_20b.txt', 
            'textfile_40b.txt', 
            'textfile_80b.txt', 
            'textfile_160b.txt', 
            'textfile_320b.txt', 
            'textfile_640b.txt', 
            'textfile_1280b.txt', 
            'textfile_2560b.txt', 
            'textfile_5120b.txt', 
            'textfile_10240b.txt', 
            'textfile_20480b.txt', 
            'textfile_40960b.txt'],
        assembly_structure='synthesis',
        encoding_scheme='max_density_encoding',
        library_name='',
        mean=20,
        vol=1000000 / Avogadro,
        std_dev=1,
        hybridisation_steps=10000,
        inner_error_correction='',
        outer_error_correction='',
        dna_barcode_length=34,  # in bp here not in oligos
        codeword_maxlength_positions=18,  # in bp here not in oligos
        years=0,
        storage_conditions=None,
        codeword_length=501,  # in bp here not in oligos
        percent_of_symbols=2,
        index_carry_length=34,  # in bp here not in oligos
        synthesis_method=68,
        sequencing_method=None,
        reed_solo_percentage=0.8,
        sigma_amount=None,
        theory='no'
    )
]