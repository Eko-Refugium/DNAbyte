"""
Example: Using encoding analysis to tune simulation parameters

This demonstrates the workflow:
1. Analyze each encoding to get metrics
2. Use metrics to set appropriate simulation parameters
3. Run benchmarks with tuned configurations
"""

from dnabyte.params import Params
from simulations.encoding_analysis import AnalyzeEncoding


def analyze_and_tune_encoding(encoding_method, base_params):
    """
    Analyze an encoding and return tuned parameters for simulation.
    
    :param encoding_method: Name of encoding (e.g., 'max_density', 'church')
    :param base_params: Base Params object with common settings
    :return: Tuned Params object and analysis metrics
    """
    
    # Create analyzer params for this encoding
    analyzer_params = Params(
        name=f'analyze_{encoding_method}',
        filename='textfile_40b.txt',
        encoding_method=encoding_method,
        binarization_method='default',
        library=base_params.library if hasattr(base_params, 'library') else None,
    )
    
    # Run analysis
    analyzer = AnalyzeEncoding(analyzer_params)
    metrics = analyzer.calculate_all_metrics()
    
    # Print report
    analyzer.print_report()
    
    # Get recommendations
    recommendations = analyzer.recommend_parameters()
    
    print(f"Recommended parameters for {encoding_method}:")
    for key, value in recommendations.items():
        print(f"  {key}: {value}")
    print()
    
    return metrics, recommendations


def setup_tuned_simulation_params(encoding_methods, base_config):
    """
    Set up simulation parameters with encoding-specific tuning.
    
    :param encoding_methods: List of encodings to analyze and simulate
    :param base_config: Dictionary with common simulation settings
    :return: List of Params objects for simulation
    """
    
    params_list = []
    metrics_by_encoding = {}
    
    print("="*70)
    print("PHASE 1: ENCODING ANALYSIS")
    print("="*70)
    
    # Analyze each encoding and get recommendations
    for encoding_method in encoding_methods:
        print(f"\nAnalyzing {encoding_method}...")
        
        base_params = Params(
            name=f'base_{encoding_method}',
            filename='textfile_40b.txt',
            encoding_method=encoding_method,
            **base_config
        )
        
        metrics, recommendations = analyze_and_tune_encoding(
            encoding_method, base_params
        )
        metrics_by_encoding[encoding_method] = {
            'metrics': metrics,
            'recommendations': recommendations
        }
    
    print("\n" + "="*70)
    print("PHASE 2: PARAMETER TUNING")
    print("="*70)
    
    # Create tuned parameter sets based on analysis
    error_rates = [0.05, 0.10, 0.15]
    repeats = 3
    
    for encoding_method in encoding_methods:
        rec = metrics_by_encoding[encoding_method]['recommendations']
        
        print(f"\nCreating parameters for {encoding_method}")
        print(f"  Assembly probability: {rec['assembly_probability']:.2f}")
        print(f"  Mean copies: {rec['mean_copies']}")
        
        for error_rate in error_rates:
            # Skip if error rate exceeds encoding capability
            max_error = metrics_by_encoding[encoding_method]['metrics'].get('max_error_rate', 0.10)
            if error_rate > max_error:
                print(f"  Skipping {error_rate:.1%} (exceeds capability of {max_error:.1%})")
                continue
            
            for repeat in range(repeats):
                params = Params(
                    name=f'{encoding_method}_err{error_rate}_{repeat+1}',
                    filename='textfile_40b.txt',
                    encoding_method=encoding_method,
                    binarization_method='default',
                    synthesis_method=None,
                    storage_conditions='biogene',
                    sequencing_method='iid',
                    iid_substitution_rate=error_rate * 0.7,
                    iid_insertion_rate=error_rate * 0.15,
                    iid_deletion_rate=error_rate * 0.15,
                    # Use recommendations from analysis
                    assembly_probability=rec['assembly_probability'],
                    mean=rec['mean_copies'],
                    std_dev=max(1, int(rec['mean_copies'] * 0.1)),
                    years=10,
                    codeword_length=500,
                    dna_barcode_length=75,
                )
                params_list.append(params)
    
    print(f"\n{len(params_list)} parameter configurations created")
    return params_list, metrics_by_encoding


# Example usage
if __name__ == '__main__':
    encoding_methods = ['max_density', 'church', 'wukong']
    
    base_config = {
        'binarization_method': 'default',
        'synthesis_method': None,
        'storage_conditions': 'biogene',
        'sequencing_method': 'iid',
        'years': 10,
        'codeword_length': 500,
        'dna_barcode_length': 75,
    }
    
    # Run analysis and parameter tuning
    params_list, metrics = setup_tuned_simulation_params(encoding_methods, base_config)
    
    print("\n" + "="*70)
    print("Ready to run simulations with tuned parameters!")
    print("="*70)
