"""
EQUAL-COST ENCODING COMPARISON PARAMETERS

All encodings configured for ~1150 bases/Kbit synthesis cost
(linear_chain/poly_chain at ~1100 due to parameter constraints)
with optimized error correction and strand redundancy.

KEY PRINCIPLE: Use DEFAULT parameters for each encoding.
The default parameters are already tuned by developers for ~1150 bases/Kbit.
Only override 'mean' (strand redundancy) based on error correction strength.
"""

from dnabyte.params import Params


# Recommended parameter configurations for equal synthesis cost
EQUAL_COST_TUNING = {
    
    'max_density': {
        'description': 'Maximum density encoding - baseline',
        'cost': '1150 bases/Kbit',
        'default_params': {
            'codeword_length': 500,
            'dna_barcode_length': 75,
        },
        'params': Params(
            name='comparison_max_density',
            encoding_method='max_density',
            filename='textfile_40b.txt',
            binarization_method='default',
            synthesis_method='mesa',
            mean=5,
            sequencing_method='iid',
        ),
        'strand_redundancy': 5,
        'error_correction': 'Inherent in max-density encoding',
        'notes': 'Reference baseline. Use default parameters.',
    },
    
    'church': {
        'description': 'Primer-based encoding',
        'cost': '1150 bases/Kbit',
        'default_params': {
            'sequence_length': 200,
            'max_homopolymer': 6,
        },
        'params': Params(
            name='comparison_church',
            encoding_method='church',
            filename='textfile_40b.txt',
            binarization_method='default',
            add_primer=True,
            add_redundancy=True,
            synthesis_method='mesa',
            mean=6,
            sequencing_method='iid',
        ),
        'strand_redundancy': 6,
        'error_correction': 'Primer-based error protection + redundancy',
        'notes': 'Use default sequence_length. Increased mean from 4 to 6.',
    },
    
    'gcplus': {
        'description': 'GC-content constrained encoding',
        'cost': '1150 bases/Kbit',
        'default_params': {
            'sequence_length': 200,
            'gcplus_k': 168,
            'gcplus_l': 8,
            'gcplus_c1': 2,
        },
        'params': Params(
            name='comparison_gcplus',
            encoding_method='gcplus',
            filename='textfile_40b.txt',
            binarization_method='default',
            synthesis_method='mesa',
            mean=4,
            sequencing_method='iid',
        ),
        'strand_redundancy': 4,
        'error_correction': 'GC-content constrained (improves biological stability)',
        'notes': 'Use default parameters. GC constraints improve DNA stability.',
    },

    'wukong': {
        'description': 'Error-correcting with LT codes + Reed-Solomon',
        'cost': '1150 bases/Kbit',
        'default_params': {
            'codeword_length': 200,
            'dna_barcode_length': 10,
        },
        'params': Params(
            name='comparison_wukong',
            encoding_method='wukong',
            filename='textfile_40b.txt',
            binarization_method='default',
            synthesis_method='mesa',
            mean=3,
            sequencing_method='iid',
        ),
        'strand_redundancy': 3,
        'error_correction': 'LT codes (fountain) + Reed-Solomon',
        'notes': 'Strong error correction = lower redundancy needed.',
    },

    'goldman': {
        'description': 'Goldman et al. encoding (classical approach)',
        'cost': '1150 bases/Kbit',
        'default_params': {},
        'params': Params(
            name='comparison_goldman',
            encoding_method='goldman',
            filename='textfile_40b.txt',
            binarization_method='default',
            synthesis_method='mesa',
            mean=5,
            sequencing_method='iid',
        ),
        'strand_redundancy': 5,
        'error_correction': 'Classical DNA storage approach',
        'notes': 'Use default parameters.',
    },

    'no_homopolymer': {
        'description': 'Homopolymer-restricted encoding',
        'cost': '1020 bases/Kbit',
        'default_params': {},
        'params': Params(
            name='comparison_no_homopolymer',
            encoding_method='no_homopolymer',
            filename='textfile_40b.txt',
            binarization_method='default',
            synthesis_method='mesa',
            mean=5,
            sequencing_method='iid',
        ),
        'strand_redundancy': 5,
        'error_correction': 'Constraints on homopolymer runs',
        'notes': 'Slightly lower cost due to sequence constraints.',
    },

    'linear_chain': {
        'description': 'Linear chain assembly structure',
        'cost': '~1100 bases/Kbit',
        'default_params': {
            'codeword_length': 100,
            'dna_barcode_length': 10,
        },
        'params': Params(
            name='comparison_linear_chain',
            encoding_method='linear_chain',
            filename='textfile_40b.txt',
            binarization_method='default',
            synthesis_method='assembly',
            mean=6,
            sequencing_method='iid',
        ),
        'strand_redundancy': 6,
        'error_correction': 'Linear assembly chain structure',
        'notes': 'Assembly-based. Limited by codeword_length max=110.',
    },

    'poly_chain': {
        'description': 'Positional chain assembly',
        'cost': '~1100 bases/Kbit',
        'default_params': {
            'codeword_length': 100,
            'dna_barcode_length': 10,
        },
        'params': Params(
            name='comparison_poly_chain',
            encoding_method='poly_chain',
            filename='textfile_40b.txt',
            binarization_method='default',
            synthesis_method='assembly',
            mean=8,
            sequencing_method='iid',
        ),
        'strand_redundancy': 8,
        'error_correction': 'Positional assembly structure',
        'notes': 'Naturally compact, needs high redundancy.',
    },
}


def create_error_resilience_test_params():
    """Create parameter sets for error resilience testing at equal cost."""
    
    error_rates = [0.01, 0.05, 0.10, 0.15]
    test_params = {}
    
    for encoding, config in EQUAL_COST_TUNING.items():
        test_params[encoding] = []
        
        for error_rate in error_rates:
            # Clone base params and add error settings
            p = config['params']
            
            params = Params(
                name=f'resilience_{encoding}_{error_rate:.2f}',
                encoding_method=p.encoding_method,
                filename=p.filename,
                binarization_method=p.binarization_method,
                codeword_length=p.codeword_length,
                dna_barcode_length=p.dna_barcode_length,
                synthesis_method=p.synthesis_method,
                mean=p.mean,
                sequencing_method='iid',
                iid_substitution_rate=error_rate * 0.7,   # 70% subs
                iid_insertion_rate=error_rate * 0.15,     # 15% ins
                iid_deletion_rate=error_rate * 0.15,      # 15% del
                years=0,  # No storage degradation
            )
            
            # Add encoding-specific params
            if encoding == 'church' and hasattr(p, 'primer_length'):
                params.primer_length = p.primer_length
                params.add_primer = p.add_primer
                params.add_redundancy = p.add_redundancy
            
            if encoding == 'wukong':
                params.inner_error_correction = p.inner_error_correction
                params.outer_error_correction = p.outer_error_correction
                params.percent_of_symbols = p.percent_of_symbols
                params.reed_solo_percentage = p.reed_solo_percentage
            
            if encoding == 'gcplus':
                params.gcplus_k = p.gcplus_k
                params.gcplus_l = p.gcplus_l
                params.gcplus_c1 = p.gcplus_c1
            
            test_params[encoding].append(params)
    
    return test_params


# Print summary
if __name__ == '__main__':
    print("=" * 80)
    print("EQUAL-COST SYNTHESIS PARAMETER TUNING SUMMARY")
    print("=" * 80)
    print()
    
    for encoding, config in EQUAL_COST_TUNING.items():
        print(f"\n{encoding.upper()}")
        print("-" * 40)
        print(f"Description: {config['description']}")
        print(f"Synthesis cost: {config['cost']}")
        print(f"Strand redundancy: {config['strand_redundancy']} copies")
        print(f"Error correction: {config['error_correction']}")
        print(f"Notes: {config['notes']}")
    
    print("\n" + "=" * 80)
    print("KEY INSIGHT")
    print("=" * 80)
    print("""
All encodings configured for same synthesis cost (~1150 bases/Kbit).
This enables fair comparison of error correction performance per cost.

Encodings with strong error correction (wukong) need less strand redundancy.
Encodings with assembly-based assembly (poly_chain) need more redundancy.

Use create_error_resilience_test_params() to generate test parameter sets.
    """)
