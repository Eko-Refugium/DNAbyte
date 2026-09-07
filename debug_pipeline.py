"""
Debug script for DNA encoding pipeline
Traces execution step-by-step to identify failure points in the recovery process
"""
import os
import sys
import traceback
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

from dnabyte.data_classes.insilicodna import InSilicoDNA
from dnabyte.encode import Encode
from dnabyte.params import Params
from dnabyte.sequence import SimulateSequencing
from dnabyte.store import SimulateStorage
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.binarize import Binarize
from dnabyte.data_classes import Data
from dnabyte.cluster import Cluster
from dnabyte.consensus import Consensus


def print_section(title):
    """Print a formatted section header"""
    print("\n" + "="*80)
    print(f"  {title}")
    print("="*80)


def print_step(step_num, description):
    """Print a step marker"""
    print(f"\n[STEP {step_num}] {description}")
    print("-" * 60)


def debug_encoding_pipeline(encoding, params, file_path='./tests/testfiles/testfilesimsall.txt', run_num=1):
    """
    Debug a single encoding through the entire pipeline
    
    Args:
        encoding: Name of encoding method
        params: Params object for the encoding
        file_path: Path to test file
        run_num: Run number for identification
    """
    print_section(f"DEBUGGING {encoding.upper()} - Run {run_num}")
    
    try:
        # STEP 1: Binarization
        print_step(1, "BINARIZATION")
        print(f"  File path: {file_path}")
        print(f"  Binarization method: {params.binarization_method}")
        
        binarizer = Binarize(params)
        data_obj = Data(file_paths=[file_path])
        binary_code = binarizer.binarize(data_obj)
        
        print(f"  ✓ Binary code length: {len(binary_code.data) if hasattr(binary_code, 'data') else 'N/A'}")
        print(f"  ✓ Binary code type: {type(binary_code)}")
        if hasattr(binary_code, 'data') and binary_code.data:
            print(f"  ✓ Sample (first 50 bits): {str(binary_code.data)[:50]}...")
        
        # STEP 2: Encoding
        print_step(2, "ENCODING")
        print(f"  Encoding method: {encoding}")
        print(f"  Encoding params: {vars(params)}")
        
        coder = Encode(params)
        data_enc, info = coder.encode(binary_code)
        
        print(f"  ✓ Encoded sequences count: {len(data_enc.data)}")
        print(f"  ✓ Encoded data type: {type(data_enc)}")
        if data_enc.data:
            print(f"  ✓ Sample sequence lengths: {[len(seq) for seq in data_enc.data[:5]]}")
            print(f"  ✓ Total encoded bp: {sum(len(seq) for seq in data_enc.data)}")
        print(f"  ℹ Info: {info}")
        
        # STEP 3: Synthesis Simulation
        print_step(3, "SYNTHESIS SIMULATION")
        print(f"  Synthesis method: {params.synthesis_method}")
        
        syn = SimulateSynthesis(params)
        data_syn, info = syn.simulate(data_enc)
        
        print(f"  ✓ Synthesized sequences count: {len(data_syn.data)}")
        print(f"  ✓ Synthesized data type: {type(data_syn)}")
        if data_syn.data:
            print(f"  ✓ Sample sequence lengths: {[len(seq) for seq in data_syn.data[:5]]}")
        print(f"  ℹ Info: {info}")
        
        # STEP 4: Storage Simulation
        print_step(4, "STORAGE SIMULATION")
        print(f"  Storage conditions: {params.storage_conditions}")
        
        sto = SimulateStorage(params)
        data_sto, info = sto.simulate(data_syn)
        
        print(f"  ✓ Stored sequences count: {len(data_sto.data)}")
        print(f"  ✓ Stored data type: {type(data_sto)}")
        if data_sto.data:
            print(f"  ✓ Sample sequence lengths: {[len(seq) for seq in data_sto.data[:5]]}")
        print(f"  ℹ Info: {info}")
        
        # STEP 5: Sequencing Simulation
        print_step(5, "SEQUENCING SIMULATION")
        print(f"  Sequencing method: {params.sequencing_method}")
        print(f"  Error rate: {params.iid_error_rate}")
        print(f"  Substitution rate: {params.iid_substitution_rate}")
        print(f"  Insertion rate: {params.iid_insertion_rate}")
        print(f"  Deletion rate: {params.iid_deletion_rate}")
        
        seq = SimulateSequencing(params)
        data_seq, info = seq.simulate(data_sto)
        
        print(f"  ✓ Sequenced reads count: {len(data_seq.data)}")
        print(f"  ✓ Sequenced data type: {type(data_seq)}")
        if data_seq.data:
            print(f"  ✓ Sample read lengths: {[len(read) for read in data_seq.data[:5]]}")
            print(f"  ✓ Total sequenced reads: {len(data_seq.data)}")
        print(f"  ℹ Info: {info}")
        
        # STEP 6: Recovery/Clustering & Consensus (if applicable)
        print_step(6, "RECOVERY PROCESS")
        print(f"  Clustering method: {params.clustering_method}")
        print(f"  Recovery method: {params.recovery_method}")
        print(f"  Expected clusters: {len(data_syn.data)} (original synthesized strands)")
        print(f"  Actual sequenced reads: {len(data_seq.data)}")
        
        if params.clustering_method is not None and params.recovery_method is not None:
            print("  Using clustering + consensus approach...")
            try:
                cluster_obj = Cluster(params)
                print(f"  → Clustering...")
                data_cluster, info = cluster_obj.cluster(data_seq)
                num_clusters = len(data_cluster.data) if hasattr(data_cluster, 'data') else 'N/A'
                print(f"    ✓ Clustered groups: {num_clusters}")
                
                # Check if clustering produced wrong number of clusters
                if isinstance(num_clusters, int) and num_clusters != len(data_syn.data):
                    print(f"    ⚠️  WARNING: Expected {len(data_syn.data)} clusters but got {num_clusters}")
                    print(f"    Ratio: {num_clusters / len(data_syn.data):.2f}x")
                
                print(f"    ℹ Clustering info: {info}")
                
                consensus_obj = Consensus(params)
                print(f"  → Calling consensus...")
                data_cor, info = consensus_obj.call(data_cluster)
                print(f"    ✓ Consensus sequences: {len(data_cor.data)}")
                print(f"    ℹ Consensus info: {info}")
            except Exception as e:
                print(f"    ✗ Clustering/Consensus failed: {type(e).__name__}: {e}")
                traceback.print_exc()
                raise
        else:
            print("  Using encoding-specific process...")
            try:
                data_cor, info = coder.process(data_seq)
                print(f"  ✓ Processed sequences: {len(data_cor.data)}")
                print(f"  ℹ Process info: {info}")
            except Exception as e:
                print(f"  ✗ Encoding-specific process failed: {type(e).__name__}: {e}")
                traceback.print_exc()
                raise
        
        # STEP 7: Decoding
        print_step(7, "DECODING")
        print(f"  Decoder type: {type(coder)}")
        
        data_dec, valid, info = coder.decode(data_cor)
        
        print(f"  ✓ Decoded data type: {type(data_dec)}")
        print(f"  ✓ Decoded valid: {valid}")
        if not valid:
            print(f"  ⚠ WARNING: Decoded data marked as invalid!")
        print(f"  ℹ Decode info: {info}")
        
        # STEP 8: Comparison
        print_step(8, "COMPARISON WITH ORIGINAL")
        
        comparison, res = data_dec.compare(data_dec, binary_code)
        
        print(f"  Comparison result: {comparison}")
        print(f"  ✓ SUCCESS!" if comparison == 'SUCCESS' else f"  ✗ FAILED: {comparison}")
        if res:
            print(f"  ℹ Comparison details: {len(data_dec.data), len(binary_code.data)}")
        
        print_section(f"✓ {encoding.upper()} - PASSED")
        return True, None
        
    except Exception as e:
        print_section(f"✗ {encoding.upper()} - FAILED")
        print(f"Error: {type(e).__name__}: {e}")
        print("\nFull traceback:")
        traceback.print_exc()
        return False, str(e)


def create_debug_params(encoding, base_params):
    """Create Params object for a specific encoding"""
    def default_values(enc):
        if enc == 'yinyang':
            return {
                'name': 'yy_default',
                'sequence_length': 200,
                'encoding_method': 'yinyang',
                'max_homopolymer': 3,
                'max_content': 0.6,
            }
        elif enc == 'goldman':
            return {
                'name': 'goldman_default',
                'sequence_length': 200,
                'encoding_method': 'goldman',
                'kmer_size_cluster': 180,
            }
        elif enc == 'church':
            return {
                'name': 'church_default',
                'sequence_length': 200,
                'encoding_method': 'church',
                'rs_num': 12,
                'max_homopolymer': 3,
                'add_redundancy': True,
                'add_primer': False,
            }
        elif enc == 'gcplus':
            return {
                'name': 'gcplus_default',
                'sequence_length': 200,
                'encoding_method': 'gcplus',
                'gcplus_k': 168,
                'gcplus_l': 8,
                'gcplus_c1': 20,
                'barcode_length': 0,
                'left_primer': '',
                'right_primer': '',
                'kmer_size_debruijn': 90,
                'kmer_size_cluster': 20,
                'kmer_threshold': 0.9,
            }
        elif enc == 'max_density':
            return {
                'name': 'max_density_default',
                'codeword_length': 200,
                'sequence_length': 200,
                'encoding_method': 'max_density',
                'max_homopolymer': 3,
                'dna_barcode_length': 8,
                'codeword_maxlength_positions': 8,
                'outer_error_correction': 'reed_solomon',
                'inner_error_correction': None,
                'reed_solo_percentage': 0.5,
                'percent_of_symbols': 1,
                'ltcode_header': 20,
            }
        elif enc == 'no_homopolymer':
            return {
                'name': 'no_homopolymer_default',
                'codeword_length': 200,
                'sequence_length': 200,
                'encoding_method': 'no_homopolymer',
                'max_homopolymer': 3,
                'dna_barcode_length': 8,
                'codeword_maxlength_positions': 8,
                'outer_error_correction': 'reed_solomon',
                'inner_error_correction': None,
                'reed_solo_percentage': 0.5,
                'percent_of_symbols': 1,
                'ltcode_header': 20,
            }
        elif enc == 'wukong':
            return {
                'name': 'wukong_default',
                'sequence_length': 200,
                'encoding_method': 'wukong',
                'max_homopolymer': 3,
                'min_gc_content': 0.4,
                'max_gc_content': 0.6,
                'rule_num': 1,
                'rs_num': 12,
                'add_redundancy': True,
                'add_primer': False,
            }
        elif enc == 'hedges':
            return {
                'name': 'hedges_default',
                'sequence_length': 200,
                'encoding_method': 'hedges',
                'max_homopolymer': 3,
                'gc_max': 8,
                'gc_window': 12,
                'hedges_coderate': 6,
                'kmer_size_cluster': 180,
            }
    
    params = base_params.copy()
    params.update(default_values(encoding))
    return Params(**params)


if __name__ == '__main__':
    import sys
    
    # Base parameters
    base_params = {
        'filename': 'testfilesimsall.txt',
        'binarization_method': 'default',
        'synthesis_method': 'nosynthpoly',
        'sequencing_method': 'iid',
        'recovery_method': 'simple',
        'min_coverage': 1,
        'clustering_method': 'similarity_cluster',
        'storage_conditions': None,
        'kmer_seed': 42,
        'mean': 1,
        'std_dev': 0,
        "iid_error_rate": 0.0,
        "iid_substitution_rate": 0.0,
        "iid_insertion_rate": 0.0,
        "iid_deletion_rate": 0.0,
        # Similarity clustering parameters (for similarity_cluster method)
        'similarity_threshold': 0.75,  # 75% similarity for 10% error rate - can override with arg 5
        'gap_penalty': 1.0,
    }
    
    # Get encoding to debug from command line or use default
    encoding_to_debug = sys.argv[1] if len(sys.argv) > 1 else 'goldman'
    error_rate = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    runs = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    clustering_method = sys.argv[4] if len(sys.argv) > 4 else 'similarity_cluster'
    similarity_threshold = float(sys.argv[5]) if len(sys.argv) > 5 else 1
    
    print(f"\n🔍 DEBUG MODE: Encoding={encoding_to_debug}, Error Rate={error_rate}, Runs={runs}")
    print(f"   Clustering Method: {clustering_method}")
    
    # Override clustering method if specified
    base_params['clustering_method'] = clustering_method
    
    # Set error rates in params
    base_params['iid_error_rate'] = error_rate
    base_params['iid_substitution_rate'] = error_rate
    base_params['iid_insertion_rate'] = 0.0
    
    # Adjust similarity threshold based on error rate (same formula as Cost_analysis_myway.py)
    if similarity_threshold is None:
        # Dynamic adjustment: 1.0 at 0% error, lower threshold at higher error rates
        similarity_threshold = max(0.85, 1.0 - (error_rate * 1.5))
    
    base_params['similarity_threshold'] = similarity_threshold
    print(f"   Similarity Threshold: {similarity_threshold:.3f} (for error_rate={error_rate})")
    print("ℹ️  Using parameters from Cost_analysis_myway.py\n")
    base_params['iid_deletion_rate'] = 0.0
    
    # Run debug for specified encoding and number of runs
    params = create_debug_params(encoding_to_debug, base_params)
    
    for run in range(1, runs + 1):
        success, error = debug_encoding_pipeline(
            encoding_to_debug, 
            params, 
            file_path='./tests/testfiles/testfilesimsall.txt',
            run_num=run
        )
        
        if not success and run < runs:
            print(f"\n⚠ Run {run} failed. Continuing to next run...\n")
    
    print("\n" + "="*80)
    print("  DEBUG SESSION COMPLETE")
    print("="*80 + "\n")
