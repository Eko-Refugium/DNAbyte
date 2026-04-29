import unittest
import os
import random
from dnabyte.params import Params
from tests.testbase_end2end_newdata import TestBase

BINARIZATION_METHODS = ['default', 'text']
STORAGE_CONDITIONS = ['biogene', 'newstorage', 'permafrost', 'random', 'roomtemperature']
SEQUENCING_METHODS = ['iid', 'illumina', 'kmere', 'mesa', 'nanopore']
SYNTHESIS_METHODS = ['mesa', 'nosynthpoly']
ERROR_METHODS = ['homopolymer_issue', 'iid']

ENCODING_METHODS = [
    'church',
    'dna_aeon',
    'gcplus',
    'goldman',
    'linear_binom',
    'linear_chain',
    'max_density',
    'no_homopolymer',
    'poly_binom',
    'poly_chain',
    'wukong',
]

line_binom_library = '20bp_Lib.csv'
poly_binom_library = 'lib_positional_22(173)m_35(2)g_24(39)p.csv'

TEST_FILE = os.path.join(os.path.dirname(__file__), '../testfiles/textfile_40b.txt')


def random_subparams(encoding):
    params = {}
    # Add assembly/percentage/symbols parameters for relevant encodings
    if encoding in ['linear_binom', 'linear_chain', 'poly_binom', 'poly_chain']:
        # Assembly-specific parameters
        params['assembly_structure'] = random.choice(['linear_assembly', 'synthesis', 'mesa'])
        params['codeword_length'] = random.choice([60, 80, 100, 120])
        params['dna_barcode_length'] = random.choice([2, 4])
        params['codeword_maxlength_positions'] = random.choice([2, 4])
        params['sigma_amount'] = random.choice([1, 2, 3])
        params['mean'] = random.choice([5, 10, 20])
        params['std_dev'] = random.choice([0])
        # Add percent and symbols for poly_binom/poly_chain
 
    # Choose safe parameters for sequence_length, primer_length, rs_num, redundancy
    params['mean']= 20
    safe_lengths = [l for l in [60, 80, 100, 120, 140] if l > 2*10]  # min primer_length is 10
    params['sequence_length'] = random.choice(safe_lengths)
    params['max_homopolymer'] = random.choice([3, 4, 5])
    params['rs_num'] = random.choice([0, 1, 2])
    params['add_redundancy'] = random.choice([False, True])
    params['years'] = random.choice([0, 1, 2, 5, 10, 20])

    params['iid_error_rate'] = round(random.uniform(0.0001, 0.01), 4)
    params['illumina_error_rate'] = round(random.uniform(0.0001, 0.01), 4)
    params['nanopore_error_rate'] = round(random.uniform(0.0001, 0.02), 4)
    params['homopolymer_error_rate'] = round(random.uniform(0.0001, 0.01), 4)

    params['text_encoding'] = random.choice(['utf-8', 'ascii'])
    params['assembly_structure'] = random.choice(['synthesis', 'mesa', 'none'])

    if params['assembly_structure'] in ['synthesis', 'mesa']:
        params['sequence_length'] = random.choice([20, 25, 30, 35, 40])

    if encoding in ['wukong', 'gcplus', 'max_density', 'poly_binom', 'poly_chain']:
        params['min_gc'] = round(random.uniform(0.3, 0.4), 2)
        params['max_gc'] = round(random.uniform(0.45, 0.55), 2)

    if encoding == 'wukong':
        params['rule_num'] = random.choice([1, 2, 3])
        params['add_primer'] = random.choice([False, True])
        # Ensure primer_length is valid and fits sequence_length
        max_primer = min(16, (params['sequence_length']//2)-1)
        primer_choices = [l for l in [10, 12, 14, 16] if l <= max_primer]
        params['primer_length'] = random.choice(primer_choices) if primer_choices else 10
        # Enforce sequence_length > 2 * primer_length if add_primer is True
        if params['add_primer']:
            min_seq_len = 2 * params['primer_length'] + 1
            if params['sequence_length'] <= 2 * params['primer_length']:
                # Pick a safe sequence_length
                params['sequence_length'] = min_seq_len + random.choice([0, 2, 4, 6, 8, 10])

    if encoding == 'dna_aeon':
        # Enforce a safe minimum chunk size (e.g. 16)
        params['dna_aeon_chunk_size'] = random.choice([16, 20, 24, 32])
        params['dna_aeon_overhead'] = round(random.uniform(0.3, 0.5), 2)
        params['dna_aeon_insert_header'] = random.choice([False, True])
        params['dna_aeon_error_correction'] = random.choice(['crc', 'nocode', 'reedsolomon', 'dna_reedsolomon'])
        params['dna_aeon_repair_symbols'] = random.choice([1, 2, 3])
        params['dna_aeon_use_dna_rules'] = random.choice([True, False])
        params['dna_aeon_drop_upper_bound'] = round(random.uniform(0.3, 0.7), 2)

    if encoding == 'gcplus':
        params['gcplus_k'] = random.choice([64, 128, 168])
        params['gcplus_l'] = random.choice([4, 8, 12])
        params['gcplus_c1'] = random.choice([1, 2, 3])

    if encoding == 'max_density':
        params['codeword_length'] = random.choice([200, 300, 400])
        params['dna_barcode_length'] = random.choice([15, 20, 25])
        params['index_carry_length'] = random.choice([10, 15])
        params['ltcode_header'] = random.choice([10, 15])
        params['sigma_amount'] = random.choice([1, 2, 3])
        params['outer_error_correction'] = random.choice([None, 'reedsolomon'])
        params['inner_error_correction'] = random.choice([None, 'ltcode'])
        params['reed_solo_percentage'] = round(random.uniform(0.7, 0.95), 2)

    if encoding == 'no_homopolymer':
        params['codeword_length'] = random.choice([200, 300, 400])
        params['dna_barcode_length'] = random.choice([15, 20, 25])
        params['index_carry_length'] = random.choice([10, 15])
        params['ltcode_header'] = random.choice([10, 15])
        params['sigma_amount'] = random.choice([1, 2, 3])
        params['outer_error_correction'] = random.choice([None, 'reedsolomon'])
        params['inner_error_correction'] = random.choice([None, 'ltcode'])
        params['reed_solo_percentage'] = round(random.uniform(0.7, 0.95), 2)

    return params


def make_test_case(params, name):
    """
    Proper unittest-compatible test wrapper.
    """

    class DynamicTest(TestBase):
        def setUp(self):
            super().setUp()
            self.params = params  # injected per-instance safely

    DynamicTest.__name__ = name
    return DynamicTest('test_logic')


def load_tests(loader, tests, pattern):
    N = 24
    suite = unittest.TestSuite()

    for i in range(N):
        encoding = random.choice(ENCODING_METHODS)
        binarization = random.choice(BINARIZATION_METHODS)
        synthesis = random.choice(SYNTHESIS_METHODS)
        storage = random.choice(STORAGE_CONDITIONS)
        sequencing = random.choice(SEQUENCING_METHODS)
        error = random.choice(ERROR_METHODS)

        subparams = random_subparams(encoding)
        extra_kwargs = {}

        if encoding in ['linear_binom', 'linear_chain', 'poly_binom', 'poly_chain']:
            synthesis = 'assembly'

        if encoding in ['linear_binom', 'linear_chain']:
            extra_kwargs['library_name'] = line_binom_library

        if encoding in ['poly_binom', 'poly_chain']:
            extra_kwargs['library_name'] = poly_binom_library

        params = Params(
            encoding_method=encoding,
            binarization_method=binarization,
            synthesis_method=synthesis,
            storage_conditions=storage,
            sequencing_method=sequencing,
            error_methods=error,
            debug=False,
            filename="test.txt",
            file_paths=["test.txt"],
            **subparams,
            **extra_kwargs
        )

        params.name = f"enc={encoding}|bin={binarization}|syn={synthesis}|sto={storage}|seq={sequencing}|err={error}"

        suite.addTest(make_test_case(params, f"TestBase_{i}"))

    return suite

def load_tests_errorfre(loader, tests, pattern):
    N = 24
    suite = unittest.TestSuite()

    for i in range(N):
        encoding = random.choice(ENCODING_METHODS)
        binarization = random.choice(BINARIZATION_METHODS)
        synthesis = None
        storage = None
        sequencing = None
        error = None

        subparams = random_subparams(encoding)
        extra_kwargs = {}

        if encoding in ['linear_binom', 'linear_chain', 'poly_binom', 'poly_chain']:
            synthesis = 'assembly'

        if encoding in ['linear_binom', 'linear_chain']:
            extra_kwargs['library_name'] = line_binom_library

        if encoding in ['poly_binom', 'poly_chain']:
            extra_kwargs['library_name'] = poly_binom_library

        params = Params(
            encoding_method=encoding,
            binarization_method=binarization,
            synthesis_method=synthesis,
            storage_conditions=storage,
            sequencing_method=sequencing,
            error_methods=error,
            debug=False,
            filename="test.txt",
            file_paths=["test.txt"],
            **subparams,
            **extra_kwargs
        )

        params.name = f"enc={encoding}|bin={binarization}|syn={synthesis}|sto={storage}|seq={sequencing}|err={error}"

        suite.addTest(make_test_case(params, f"TestBase_{i}"))

    return suite

if __name__ == '__main__':
    print("\n[INFO] Running all parameterized end-to-end tests:\n")

    suite = load_tests_errorfre(None, None, None)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    # --- Summary report ---
    total = result.testsRun
    failures = len(result.failures)
    errors = len(result.errors)
    passed = total - failures - errors

    print("\n================ TEST SUMMARY ================\n")
    print(f"Total tests run : {total}")
    print(f"Passed          : {passed}")
    print(f"Failures        : {failures}")
    print(f"Errors          : {errors}")
    print("\n==============================================\n")

    import sys
    sys.exit(not result.wasSuccessful())