"""
Error Summary Generator

Scans all test log files and extracts errors into a single summary file.
Run this after tests complete to get a consolidated error report.
"""

import os
import re
from datetime import datetime
from pathlib import Path


def extract_errors_from_logs(log_dir='tests/testlogs'):
    """
    Scans all log files in the testlogs directory and extracts errors.
    Returns a dictionary organized by error type and test.
    """
    errors = {
        'decoding_errors': [],
        'encoding_errors': [],
        'parameter_errors': [],
        'other_errors': [],
        'failed_tests': [],
        'file_recreation_errors': [],  # Non-critical: file restoration issues
    }
    
    if not os.path.exists(log_dir):
        return errors
    
    log_files = sorted([f for f in os.listdir(log_dir) if f.endswith('.log')])
    
    for log_file in log_files:
        log_path = os.path.join(log_dir, log_file)
        
        try:
            with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                lines = content.split('\n')
                
            # Extract test name
            test_name = log_file
            test_setting = None
            for line in lines:
                if 'TEST SETTING:' in line:
                    test_setting = line.split('TEST SETTING:')[-1].strip()
                    break
            
            # Look for ERROR lines
            in_step10 = False
            for i, line in enumerate(lines):
                # Track if we're in STEP10 (file restoration)
                if 'STEP10:' in line or 'RESTORE DATA' in line:
                    in_step10 = True
                elif 'STEP' in line and 'STEP10' not in line:
                    in_step10 = False
                
                if ' - ERROR - ' in line:
                    error_text = line.split(' - ERROR - ')[-1].strip()
                    
                    # File recreation errors are non-critical
                    if in_step10 or 'debinarize' in error_text.lower() or 'restore' in error_text.lower():
                        errors['file_recreation_errors'].append({
                            'test_file': test_name,
                            'test_setting': test_setting,
                            'error': error_text,
                        })
                    # Categorize other errors
                    elif 'decode' in error_text.lower() or 'decoding' in error_text.lower():
                        errors['decoding_errors'].append({
                            'test_file': test_name,
                            'test_setting': test_setting,
                            'error': error_text,
                        })
                    elif 'encod' in error_text.lower():
                        errors['encoding_errors'].append({
                            'test_file': test_name,
                            'test_setting': test_setting,
                            'error': error_text,
                        })
                    elif 'param' in error_text.lower() or 'chunk' in error_text.lower():
                        errors['parameter_errors'].append({
                            'test_file': test_name,
                            'test_setting': test_setting,
                            'error': error_text,
                        })
                    else:
                        errors['other_errors'].append({
                            'test_file': test_name,
                            'test_setting': test_setting,
                            'error': error_text,
                        })
                
                # Look for test failures (but skip STEP10 failures)
                if not in_step10 and ('AssertionError:' in line or 'STATUS: ERROR' in line):
                    errors['failed_tests'].append({
                        'test_file': test_name,
                        'test_setting': test_setting,
                        'error': line.strip(),
                    })
                    
        except Exception as e:
            print(f"Failed to process {log_file}: {e}")
    
    return errors


def generate_summary_file(output_file='tests/ERROR_SUMMARY.txt'):
    """
    Generates a consolidated error summary file grouped by parameters.
    """
    errors = extract_errors_from_logs()
    
    # Reorganize all errors by test_setting (parameters)
    errors_by_params = {}
    
    for category in errors:
        for error_info in errors[category]:
            test_setting = error_info.get('test_setting', 'Unknown')
            if test_setting not in errors_by_params:
                errors_by_params[test_setting] = {
                    'decoding': [],
                    'encoding': [],
                    'parameter': [],
                    'other': [],
                    'failed': [],
                    'file_recreation': [],
                }
            
            if category == 'decoding_errors':
                errors_by_params[test_setting]['decoding'].append(error_info)
            elif category == 'encoding_errors':
                errors_by_params[test_setting]['encoding'].append(error_info)
            elif category == 'parameter_errors':
                errors_by_params[test_setting]['parameter'].append(error_info)
            elif category == 'other_errors':
                errors_by_params[test_setting]['other'].append(error_info)
            elif category == 'failed_tests':
                errors_by_params[test_setting]['failed'].append(error_info)
            elif category == 'file_recreation_errors':
                errors_by_params[test_setting]['file_recreation'].append(error_info)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("TEST ERROR SUMMARY REPORT (GROUPED BY PARAMETERS)\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")
        
        # Count summary - separate critical from non-critical
        critical_errors = sum(len(v) for k, v in errors.items() 
                            if isinstance(v, list) and k != 'file_recreation_errors')
        non_critical_errors = len(errors['file_recreation_errors'])
        
        f.write(f"CRITICAL ERRORS: {critical_errors}\n")
        f.write(f"NON-CRITICAL ERRORS (File Recreation): {non_critical_errors}\n")
        f.write(f"TOTAL ISSUES: {critical_errors + non_critical_errors}\n\n")
        
        # Quick summary by encoding
        encoding_summary = {}
        for test_setting in errors_by_params:
            # Extract encoding if available
            if '|' in test_setting:
                parts = test_setting.split('|')
                for part in parts:
                    if part.startswith('enc='):
                        encoding = part.split('=')[1]
                        if encoding not in encoding_summary:
                            encoding_summary[encoding] = 0
                        param_errors = sum(len(v) for k, v in errors_by_params[test_setting].items() 
                                         if k != 'file_recreation' and isinstance(v, list))
                        encoding_summary[encoding] += param_errors
                        break
        
        if encoding_summary:
            f.write("=" * 80 + "\n")
            f.write("ERRORS BY ENCODING TYPE\n")
            f.write("=" * 80 + "\n")
            for encoding in sorted(encoding_summary.keys()):
                count = encoding_summary[encoding]
                f.write(f"{encoding}: {count} error(s)\n")
            f.write("\n")
        
        # Now show details grouped by parameters
        f.write("=" * 80 + "\n")
        f.write("DETAILED ERRORS GROUPED BY PARAMETERS\n")
        f.write("=" * 80 + "\n\n")
        
        for test_setting in sorted(errors_by_params.keys()):
            param_errors = errors_by_params[test_setting]
            total_param_errors = sum(len(v) for k, v in param_errors.items() 
                                    if isinstance(v, list) and k != 'file_recreation')
            
            if total_param_errors > 0:
                f.write("\n" + "=" * 80 + "\n")
                f.write(f"PARAMETERS: {test_setting}\n")
                f.write("=" * 80 + "\n")
                f.write(f"Total Critical Errors: {total_param_errors}\n\n")
                
                # Encoding errors
                if param_errors['encoding']:
                    f.write("\n[ENCODING ERRORS]\n")
                    f.write("-" * 80 + "\n")
                    for error_info in param_errors['encoding']:
                        f.write(f"File: {error_info['test_file']}\n")
                        f.write(f"Error: {error_info['error']}\n")
                        f.write("-" * 40 + "\n")
                
                # Decoding errors
                if param_errors['decoding']:
                    f.write("\n[DECODING ERRORS]\n")
                    f.write("-" * 80 + "\n")
                    for error_info in param_errors['decoding']:
                        f.write(f"File: {error_info['test_file']}\n")
                        f.write(f"Error: {error_info['error']}\n")
                        f.write("-" * 40 + "\n")
                
                # Parameter errors
                if param_errors['parameter']:
                    f.write("\n[PARAMETER ERRORS]\n")
                    f.write("-" * 80 + "\n")
                    for error_info in param_errors['parameter']:
                        f.write(f"File: {error_info['test_file']}\n")
                        f.write(f"Error: {error_info['error']}\n")
                        f.write("-" * 40 + "\n")
                
                # Other errors
                if param_errors['other']:
                    f.write("\n[OTHER ERRORS]\n")
                    f.write("-" * 80 + "\n")
                    for error_info in param_errors['other']:
                        f.write(f"File: {error_info['test_file']}\n")
                        f.write(f"Error: {error_info['error']}\n")
                        f.write("-" * 40 + "\n")
                
                # Failed tests
                if param_errors['failed']:
                    f.write("\n[FAILED TESTS]\n")
                    f.write("-" * 80 + "\n")
                    for error_info in param_errors['failed']:
                        f.write(f"File: {error_info['test_file']}\n")
                        f.write(f"Error: {error_info['error']}\n")
                        f.write("-" * 40 + "\n")
        
        # Non-critical errors section
        if errors['file_recreation_errors']:
            f.write("\n" + "=" * 80 + "\n")
            f.write("[NON-CRITICAL] FILE RECREATION ERRORS\n")
            f.write("=" * 80 + "\n")
            f.write("These errors occur during file restoration (STEP10) and are not critical.\n")
            f.write("The pipeline succeeded up to this point.\n\n")
            
            # Group file recreation by params too
            file_errors_by_params = {}
            for error_info in errors['file_recreation_errors']:
                test_setting = error_info.get('test_setting', 'Unknown')
                if test_setting not in file_errors_by_params:
                    file_errors_by_params[test_setting] = []
                file_errors_by_params[test_setting].append(error_info)
            
            for test_setting in sorted(file_errors_by_params.keys()):
                if file_errors_by_params[test_setting]:
                    f.write(f"\nParameters: {test_setting}\n")
                    f.write("-" * 80 + "\n")
                    for error_info in file_errors_by_params[test_setting]:
                        f.write(f"File: {error_info['test_file']}\n")
                        f.write(f"Error: {error_info['error']}\n")
                        f.write("-" * 40 + "\n")
        
        if critical_errors == 0:
            f.write("\n✓ NO CRITICAL ERRORS FOUND!\n")
            if non_critical_errors > 0:
                f.write(f"({non_critical_errors} non-critical file recreation issues)\n")
    
    print(f"Error summary written to: {output_file}")
    return output_file


if __name__ == '__main__':
    generate_summary_file()
