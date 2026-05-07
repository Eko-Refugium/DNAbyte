"""
test_end2end_no_errors.py
=========================
Systematic end-to-end tests for every encoding method with NO error channels.

Goal: verify that every supported encoding can round-trip
  encode → synthesize 1 copy → process → decode → compare
correctly across a wide range of its own parameters.
If any test fails here the bug is purely in the encode/process/decode
pipeline, not in error-channel handling.

Error channels disabled on ALL tests:
  storage_conditions = None
  sequencing_method  = None
  error_methods      = None

Synthesis:
  nosynthpoly / mean=1, std_dev=0  – synthesis-based encodings
  assembly                          – assembly-based encodings

Total tests: ~90
"""

import unittest
from dnabyte.params import Params
from tests.testbase_end2end_newdata import TestBase

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

_NO_ERR = dict(
    storage_conditions=None,
    sequencing_method=None,
    error_methods=None,
    binarization_method='default',
    filename='textfile_40b.txt',
)
_S1 = dict(synthesis_method='nosynthpoly', mean=1, std_dev=0)   # 1 clean copy

# For text binarizer: file_paths must be present at Params() construction time
_NO_ERR_TEXT = dict(
    storage_conditions=None,
    sequencing_method=None,
    error_methods=None,
    binarization_method='text',
    filename='textfile_40b.txt',
    file_paths=['textfile_40b.txt'],
    text_encoding='utf-8',
)


def _p(name, **kw):
    return Params(**{**_NO_ERR, **kw, 'name': name})


def _pt(name, **kw):
    """Like _p but uses the text binarizer (requires file_paths)."""
    return Params(**{**_NO_ERR_TEXT, **kw, 'name': name})


# ============================================================================
# CHURCH
# ============================================================================
church_params = [
    # --- primer length variants ---
    _p('church_primer20_hp2',        encoding_method='church', **_S1,
       sequence_length=200, max_homopolymer=2, rs_num=0,
       add_redundancy=True,  add_primer=True,  primer_length=20),

    _p('church_no_primer_hp4',       encoding_method='church', **_S1,
       sequence_length=200, max_homopolymer=4, rs_num=0,
       add_redundancy=False, add_primer=False, primer_length=0),

    # --- sequence length variants (seq >= 200 work best) ---
    _p('church_seq250_primer20',     encoding_method='church', **_S1,
       sequence_length=250, max_homopolymer=4, rs_num=0,
       add_redundancy=True,  add_primer=True,  primer_length=20),

    _p('church_seq300_primer20',     encoding_method='church', **_S1,
       sequence_length=300, max_homopolymer=4, rs_num=0,
       add_redundancy=True,  add_primer=True,  primer_length=20),

    # --- Reed-Solomon error correction ---
    _p('church_rs1_primer20',        encoding_method='church', **_S1,
       sequence_length=200, max_homopolymer=4, rs_num=1,
       add_redundancy=True,  add_primer=True,  primer_length=20),

    _p('church_rs2_primer20',        encoding_method='church', **_S1,
       sequence_length=250, max_homopolymer=4, rs_num=2,
       add_redundancy=True,  add_primer=True,  primer_length=20),

    # --- no redundancy ---
    _p('church_no_redundancy',       encoding_method='church', **_S1,
       sequence_length=200, max_homopolymer=4, rs_num=0,
       add_redundancy=False, add_primer=True,  primer_length=20),

    # --- strict vs loose homopolymer ---
    _p('church_hp6_seq200',          encoding_method='church', **_S1,
       sequence_length=200, max_homopolymer=6, rs_num=0,
       add_redundancy=True,  add_primer=True,  primer_length=20),
]

# ============================================================================
# WUKONG
# ============================================================================
wukong_params = [
    # --- rule_num variants ---
    _p('wukong_rule1_gc40_60',       encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=1, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    _p('wukong_rule2_gc40_60',       encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=2, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    _p('wukong_rule3_gc40_60',       encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=3, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    # --- GC window variants ---
    _p('wukong_gc30_70',             encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=6,
       min_gc=0.30, max_gc=0.70, rule_num=1, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    _p('wukong_gc45_55',             encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=4,
       min_gc=0.45, max_gc=0.55, rule_num=1, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    # --- sequence length variants (seq >= 200 work best) ---
    _p('wukong_seq250_primer20',     encoding_method='wukong', **_S1,
       sequence_length=250, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=1, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    _p('wukong_seq300_primer20',     encoding_method='wukong', **_S1,
       sequence_length=300, max_homopolymer=6,
       min_gc=0.30, max_gc=0.70, rule_num=1, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    # --- primer variants ---
    _p('wukong_no_primer',           encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=1, rs_num=0,
       add_redundancy=False, add_primer=False, primer_length=0),

    # --- RS error correction ---
    _p('wukong_rs1',                 encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=1, rs_num=1,
       add_redundancy=True,  add_primer=True, primer_length=20),

    _p('wukong_rs2',                 encoding_method='wukong', **_S1,
       sequence_length=250, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=1, rs_num=2,
       add_redundancy=True,  add_primer=True, primer_length=20),

    # --- no redundancy ---
    _p('wukong_no_redundancy',       encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=4,
       min_gc=0.40, max_gc=0.60, rule_num=2, rs_num=0,
       add_redundancy=False, add_primer=True, primer_length=20),

    # --- strict homopolymer ---
    _p('wukong_hp2',                 encoding_method='wukong', **_S1,
       sequence_length=200, max_homopolymer=2,
       min_gc=0.40, max_gc=0.60, rule_num=1, rs_num=0,
       add_redundancy=True,  add_primer=True, primer_length=20),

    # --- combined rule + GC stress ---
    _p('wukong_rule2_gc30_70_rs1',   encoding_method='wukong', **_S1,
       sequence_length=250, max_homopolymer=4,
       min_gc=0.30, max_gc=0.70, rule_num=2, rs_num=1,
       add_redundancy=True,  add_primer=True, primer_length=20),
]

# ============================================================================
# GOLDMAN
# ============================================================================
goldman_params = [
    _p('goldman_seq100',  encoding_method='goldman', **_S1,
       sequence_length=100, add_primer=False, primer_length=0),
    _p('goldman_seq150',  encoding_method='goldman', **_S1,
       sequence_length=150, add_primer=False, primer_length=0),
    _p('goldman_seq200',  encoding_method='goldman', **_S1,
       sequence_length=200, add_primer=False, primer_length=0),
    _p('goldman_seq250',  encoding_method='goldman', **_S1,
       sequence_length=250, add_primer=False, primer_length=0),
    _p('goldman_seq300',  encoding_method='goldman', **_S1,
       sequence_length=300, add_primer=False, primer_length=0),
]

# ============================================================================
# GC+
# ============================================================================
gcplus_params = [
    # k=168 variants (only with larger l values to avoid too-small codewords)
    _p('gcplus_k168_l8_c2',    encoding_method='gcplus', **_S1,
       sequence_length=200, gcplus_k=168, gcplus_l=8,  gcplus_c1=2),
    # k=128 variants
    _p('gcplus_k128_l8_c2',    encoding_method='gcplus', **_S1,
       sequence_length=200, gcplus_k=128, gcplus_l=8,  gcplus_c1=2),
    # k=64 variants
    _p('gcplus_k64_l8_c2',     encoding_method='gcplus', **_S1,
       sequence_length=200, gcplus_k=64,  gcplus_l=8,  gcplus_c1=2),
    # longer sequences
    _p('gcplus_k168_l12_seq300', encoding_method='gcplus', **_S1,
       sequence_length=300, gcplus_k=168, gcplus_l=12, gcplus_c1=3),
    _p('gcplus_k128_l12_seq250', encoding_method='gcplus', **_S1,
       sequence_length=250, gcplus_k=128, gcplus_l=12, gcplus_c1=2),
]

# ============================================================================
# DNA-AEON
# ============================================================================
# dna_aeon_params = [
#     # error correction modes
#     _p('dna_aeon_crc_chunk10',          encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=10,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     _p('dna_aeon_nocode_chunk20',        encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=20,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='nocode',         dna_aeon_repair_symbols=0,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     _p('dna_aeon_rs_chunk16',            encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=16,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='reedsolomon',    dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     _p('dna_aeon_dna_rs_chunk16',        encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=16,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='dna_reedsolomon', dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     # DNA rules off
#     _p('dna_aeon_crc_no_dna_rules',      encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=10,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=False, dna_aeon_drop_upper_bound=0.5),

#     # overhead variants
#     _p('dna_aeon_overhead30',            encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=10,  dna_aeon_overhead=0.30,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     _p('dna_aeon_overhead50',            encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=10,  dna_aeon_overhead=0.50,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     # chunk size variants
#     _p('dna_aeon_chunk24',               encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=24,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     _p('dna_aeon_chunk32',               encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=32,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     # repair symbol variants
#     _p('dna_aeon_rs_repair1',            encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=16,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='reedsolomon',    dna_aeon_repair_symbols=1,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     _p('dna_aeon_rs_repair3',            encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=16,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='reedsolomon',    dna_aeon_repair_symbols=3,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     # longer sequence
#     _p('dna_aeon_seq300_crc',            encoding_method='dna_aeon', **_S1,
#        sequence_length=300, dna_aeon_chunk_size=10,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),

#     # insert_header on
#     _p('dna_aeon_insert_header',         encoding_method='dna_aeon', **_S1,
#        sequence_length=200, dna_aeon_chunk_size=10,  dna_aeon_overhead=0.40,
#        dna_aeon_error_correction='crc',           dna_aeon_repair_symbols=2,
#        dna_aeon_insert_header=True,
#        dna_aeon_use_dna_rules=True,  dna_aeon_drop_upper_bound=0.5),
# ]

# ============================================================================
# MAX DENSITY
# ============================================================================
max_density_params = [
    # no error correction, various codeword sizes
    _p('maxdensity_cw200_bc30',          encoding_method='max_density', **_S1,
       codeword_length=200, dna_barcode_length=30,
       inner_error_correction=None, outer_error_correction=None),

    _p('maxdensity_cw300_bc45',          encoding_method='max_density', **_S1,
       codeword_length=300, dna_barcode_length=45,
       inner_error_correction=None, outer_error_correction=None),

    _p('maxdensity_cw400_bc60',          encoding_method='max_density', **_S1,
       codeword_length=400, dna_barcode_length=60,
       inner_error_correction=None, outer_error_correction=None),

    _p('maxdensity_cw500_bc75',          encoding_method='max_density', **_S1,
       codeword_length=500, dna_barcode_length=75,
       inner_error_correction=None, outer_error_correction=None),

    _p('maxdensity_cw700_bc100',         encoding_method='max_density', **_S1,
       codeword_length=700, dna_barcode_length=100,
       inner_error_correction=None, outer_error_correction=None),

    # RS outer correction – different percentages
    _p('maxdensity_cw400_rs70',          encoding_method='max_density', **_S1,
       codeword_length=400, dna_barcode_length=60,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.70),

    _p('maxdensity_cw400_rs80',          encoding_method='max_density', **_S1,
       codeword_length=400, dna_barcode_length=60,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.80),

    _p('maxdensity_cw500_rs90',          encoding_method='max_density', **_S1,
       codeword_length=500, dna_barcode_length=75,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.90),
]

# ============================================================================
# NO HOMOPOLYMER
# ============================================================================
no_homopolymer_params = [
    # no error correction, various codeword sizes
    _p('nohomopoly_cw200_bc30',          encoding_method='no_homopolymer', **_S1,
       codeword_length=200, dna_barcode_length=30, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction=None),

    _p('nohomopoly_cw300_bc45',          encoding_method='no_homopolymer', **_S1,
       codeword_length=300, dna_barcode_length=45, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction=None),

    _p('nohomopoly_cw400_bc60',          encoding_method='no_homopolymer', **_S1,
       codeword_length=400, dna_barcode_length=60, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction=None),

    _p('nohomopoly_cw500_bc75',          encoding_method='no_homopolymer', **_S1,
       codeword_length=500, dna_barcode_length=75, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction=None),

    _p('nohomopoly_cw700_bc100',         encoding_method='no_homopolymer', **_S1,
       codeword_length=700, dna_barcode_length=100, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction=None),

    # barcode size variant
    _p('nohomopoly_cw500_bc50',          encoding_method='no_homopolymer', **_S1,
       codeword_length=500, dna_barcode_length=50, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction=None),

    # RS outer correction
    _p('nohomopoly_cw400_rs80',          encoding_method='no_homopolymer', **_S1,
       codeword_length=400, dna_barcode_length=60, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.80),

    _p('nohomopoly_cw500_rs90',          encoding_method='no_homopolymer', **_S1,
       codeword_length=500, dna_barcode_length=75, codeword_maxlength_positions=10,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.90),
]

# ============================================================================
# LINEAR BINOM  (assembly-based)
# ============================================================================
linear_binom_params = [
    _p('linear_binom_mean5_bc2',         encoding_method='linear_binom',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=5,  std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction=None),

    _p('linear_binom_mean10_bc2',        encoding_method='linear_binom',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=10, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction=None),

    _p('linear_binom_mean10_bc4',        encoding_method='linear_binom',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=10, std_dev=0, dna_barcode_length=4,
       inner_error_correction=None, outer_error_correction=None),

    _p('linear_binom_mean20_bc2',        encoding_method='linear_binom',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=20, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction=None),

    # RS outer correction
    _p('linear_binom_mean10_bc2_rs80',   encoding_method='linear_binom',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=10, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.80),

    _p('linear_binom_mean20_bc2_rs90',   encoding_method='linear_binom',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=20, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.90),
]

# ============================================================================
# LINEAR CHAIN  (assembly-based)
# ============================================================================
linear_chain_params = [
    _p('linear_chain_mean5_bc2',         encoding_method='linear_chain',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=5,  std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction=None),

    _p('linear_chain_mean10_bc2',        encoding_method='linear_chain',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=10, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction=None),

    _p('linear_chain_mean10_bc4',        encoding_method='linear_chain',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=10, std_dev=0, dna_barcode_length=4,
       inner_error_correction=None, outer_error_correction=None),

    _p('linear_chain_mean20_bc2',        encoding_method='linear_chain',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=20, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction=None),

    # RS outer correction
    _p('linear_chain_mean10_bc2_rs80',   encoding_method='linear_chain',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=10, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.80),

    _p('linear_chain_mean20_bc2_rs90',   encoding_method='linear_chain',
       synthesis_method='assembly', library_name='20bp_Lib.csv',
       mean=20, std_dev=0, dna_barcode_length=2,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.90),
]

# ============================================================================
# POLY BINOM  (assembly-based)
# ============================================================================
poly_binom_params = [
    _p('poly_binom_mean1',               encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=1,  std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    _p('poly_binom_mean5',               encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=5,  std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    _p('poly_binom_mean10',              encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=10, std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    _p('poly_binom_mean20',              encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=20, std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    # RS outer correction
    _p('poly_binom_mean10_rs80',         encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=10, std_dev=0,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.80),

    _p('poly_binom_mean20_rs90',         encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=20, std_dev=0,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.90),

    # ltcode inner correction
    _p('poly_binom_mean10_lt',           encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=10, std_dev=0,
       inner_error_correction='ltcode',  outer_error_correction=None,
       index_carry_length=2, ltcode_header=2, percent_of_symbols=2),

    _p('poly_binom_mean10_lt_rs80',      encoding_method='poly_binom',
       synthesis_method='assembly',
       library_name='lib_positional_22(81)m_36(2)g_21(39)p.csv',
       mean=10, std_dev=0,
       inner_error_correction='ltcode',  outer_error_correction='reedsolomon',
       index_carry_length=2, ltcode_header=2, percent_of_symbols=2,
       reed_solo_percentage=0.80),
]

# ============================================================================
# POLY CHAIN  (assembly-based)
# ============================================================================
poly_chain_params = [
    _p('poly_chain_mean1',               encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=1,  std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    _p('poly_chain_mean5',               encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=5,  std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    _p('poly_chain_mean10',              encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=10, std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    _p('poly_chain_mean20',              encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=20, std_dev=0,
       inner_error_correction=None, outer_error_correction=None),

    # RS outer correction
    _p('poly_chain_mean10_rs80',         encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=10, std_dev=0,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.80),

    _p('poly_chain_mean20_rs90',         encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=20, std_dev=0,
       inner_error_correction=None, outer_error_correction='reedsolomon',
       reed_solo_percentage=0.90),

    # ltcode inner correction
    _p('poly_chain_mean10_lt',           encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=10, std_dev=0,
       inner_error_correction='ltcode',  outer_error_correction=None,
       index_carry_length=2, ltcode_header=2, percent_of_symbols=2),

    _p('poly_chain_mean10_lt_rs80',      encoding_method='poly_chain',
       synthesis_method='assembly', library_name='polymeraselibinparts.txt',
       mean=10, std_dev=0,
       inner_error_correction='ltcode',  outer_error_correction='reedsolomon',
       index_carry_length=2, ltcode_header=2, percent_of_symbols=2,
       reed_solo_percentage=0.80),
]

# ============================================================================
# Cross-cutting: text binarizer with representative encodings
# (uses _pt() which supplies the file_paths required by the text binarizer)
# ============================================================================
text_binarization_params = [
    _pt('church_text_binarize',      encoding_method='church', **_S1,
        sequence_length=200, max_homopolymer=4, rs_num=0,
        add_redundancy=True, add_primer=True, primer_length=20),

    _pt('church_text_binarize_rs1',  encoding_method='church', **_S1,
        sequence_length=200, max_homopolymer=4, rs_num=1,
        add_redundancy=True, add_primer=True, primer_length=20),

    _pt('wukong_text_binarize',      encoding_method='wukong', **_S1,
        sequence_length=200, max_homopolymer=4,
        min_gc=0.40, max_gc=0.60, rule_num=1, rs_num=0,
        add_redundancy=True, add_primer=True, primer_length=20),

    _pt('wukong_text_binarize_rule2', encoding_method='wukong', **_S1,
        sequence_length=200, max_homopolymer=4,
        min_gc=0.40, max_gc=0.60, rule_num=2, rs_num=0,
        add_redundancy=True, add_primer=True, primer_length=20),

    _pt('goldman_text_binarize',     encoding_method='goldman', **_S1,
        sequence_length=200, add_primer=False, primer_length=0),

    _pt('gcplus_text_binarize',      encoding_method='gcplus', **_S1,
        sequence_length=200, gcplus_k=168, gcplus_l=8, gcplus_c1=2),

    _pt('dna_aeon_text_binarize',    encoding_method='dna_aeon', **_S1,
        sequence_length=200, dna_aeon_chunk_size=10, dna_aeon_overhead=0.40,
        dna_aeon_error_correction='crc', dna_aeon_repair_symbols=2,
        dna_aeon_use_dna_rules=True,    dna_aeon_drop_upper_bound=0.5),

    _pt('max_density_text_binarize', encoding_method='max_density', **_S1,
        codeword_length=500, dna_barcode_length=75,
        inner_error_correction=None, outer_error_correction=None),

    _pt('max_density_text_binarize_rs', encoding_method='max_density', **_S1,
        codeword_length=500, dna_barcode_length=75,
        inner_error_correction=None, outer_error_correction='reedsolomon',
        reed_solo_percentage=0.80),

    _pt('nohomopoly_text_binarize',  encoding_method='no_homopolymer', **_S1,
        codeword_length=500, dna_barcode_length=75, codeword_maxlength_positions=10,
        inner_error_correction=None, outer_error_correction=None),

    _pt('nohomopoly_text_binarize_rs', encoding_method='no_homopolymer', **_S1,
        codeword_length=500, dna_barcode_length=75, codeword_maxlength_positions=10,
        inner_error_correction=None, outer_error_correction='reedsolomon',
        reed_solo_percentage=0.80),
]

# ============================================================================
# Master list
# ============================================================================
ALL_PARAMS = (
    church_params
    + wukong_params
    + goldman_params
    + gcplus_params
    # + dna_aeon_params
    + max_density_params
    + no_homopolymer_params
    + linear_binom_params
    + linear_chain_params
    + poly_binom_params
    + poly_chain_params
    + text_binarization_params
)

# ============================================================================
# Test-case factory
# ============================================================================

def _make_cls(params):
    class _T(TestBase):
        def __init__(self, methodName='test_logic'):
            super().__init__(methodName, params=params)
    _T.__name__ = f'Test_{params.name}'
    _T.__qualname__ = _T.__name__
    return _T


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for p in ALL_PARAMS:
        suite.addTest(_make_cls(p)('test_logic'))
    return suite


if __name__ == '__main__':
    import sys

    suite = load_tests(None, None, None)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    total  = result.testsRun
    passed = total - len(result.failures) - len(result.errors)
    print("\n================ NO-ERROR TEST SUMMARY ================")
    print(f"  Total  : {total}")
    print(f"  Passed : {passed}")
    print(f"  Failed : {len(result.failures)}")
    print(f"  Errors : {len(result.errors)}")
    print("=======================================================\n")

    # Write failed/error tests to file
    with open('failed_tests.txt', 'w') as f:
        if result.failures:
            f.write("===== FAILED TESTS =====\n")
            for test, traceback in result.failures:
                f.write(f"\n{test}:\n")
                f.write(traceback)
                f.write("\n" + "="*80 + "\n")
        
        if result.errors:
            f.write("\n===== ERROR TESTS =====\n")
            for test, traceback in result.errors:
                f.write(f"\n{test}:\n")
                f.write(traceback)
                f.write("\n" + "="*80 + "\n")
        
        if not result.failures and not result.errors:
            f.write("All tests passed!\n")

    sys.exit(not result.wasSuccessful())

