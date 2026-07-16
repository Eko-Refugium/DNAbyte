"""
Test Runner with Error Summary

Runs all-params tests and automatically generates an error summary file.
Usage: python run_tests_with_summary.py
"""

import os
import sys
import unittest
import subprocess
from datetime import datetime


def run_tests_and_summarize():
    """
    Runs the test suite and generates an error summary.
    """
    print("\n" + "=" * 80)
    print("RUNNING TEST SUITE WITH ERROR COLLECTION")
    print("=" * 80 + "\n")
    
    # Ensure testlogs directory exists
    os.makedirs('tests/testlogs', exist_ok=True)
    
    # Run the tests
    print(f"Starting test run at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Import and run test_end2end_allparams
    try:
        # Change to repo root if needed
        if not os.path.exists('tests/end2end'):
            print("ERROR: Cannot find tests directory. Make sure you're in the repo root.")
            sys.exit(1)
        
        # Run tests using unittest discovery
        loader = unittest.TestLoader()
        suite = loader.discover('tests/end2end', pattern='test_end2end_allparams.py')
        
        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)
        
        test_count = result.testsRun
        failures = len(result.failures)
        errors = len(result.errors)
        
        print(f"\n" + "=" * 80)
        print(f"TEST RUN COMPLETE")
        print(f"Tests run: {test_count}, Failures: {failures}, Errors: {errors}")
        print("=" * 80 + "\n")
        
    except Exception as e:
        print(f"Error running tests: {e}")
        import traceback
        traceback.print_exc()
    
    # Generate error summary
    print("Generating error summary...\n")
    try:
        from error_summary_generator import generate_summary_file
        summary_file = generate_summary_file()
        print(f"\n✓ Error summary generated: {summary_file}")
        print("\nYou can now review all errors in: tests/ERROR_SUMMARY.txt")
        
    except Exception as e:
        print(f"Error generating summary: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    run_tests_and_summarize()
