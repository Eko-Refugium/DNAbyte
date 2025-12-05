from dnabyte.synthesize import SimulateSynthesis
from dnabyte.params import Params
from dnabyte import NucleobaseCode, InSilicoDNA
import unittest

class TestMESA(unittest.TestCase):
    
    def test_simulate_synthesis_all_ids(self):
        """Test all valid mesa_synthesis_id values"""
        test_ids = [3, 4, 5, 6, 7, 68, 69, 70, 71]
        test_sequences = ['ACACTGATACAACTTTAGTTGTTTACACTGAGATCACTTTATATAGG', 'ACGGGGACATGTGTTAATATATACACACCATACAGCTCAATCATCTAGCATTCTCAGCA']
        
        for mesa_id in test_ids:
            with self.subTest(mesa_synthesis_id=mesa_id):
                # Initialize MESA synthesis simulator with current ID
                params = Params(
                    assembly_structure='synthesis', 
                    synthesis_method='mesa', 
                    mesa_synthesis_id=mesa_id, 
                    mean=100, 
                    std_dev=5)
                simulator = SimulateSynthesis(params=params)
                
                data = NucleobaseCode(test_sequences)
                
                # Simulate synthesis
                synthesized_sequences, info = simulator.simulate(data)
                
                print(synthesized_sequences)

                # Check that the output is InSilicoDNA
                self.assertIsInstance(synthesized_sequences, InSilicoDNA)
                
                # Check that the number of synthesized sequences is greater than or equal to input
                self.assertGreaterEqual(len(synthesized_sequences.data), len(data.data))
                
                # Check that info dictionary contains expected keys
                self.assertIn('average_copy_number', info)
                self.assertIn('number_of_synthesis_errors', info)
                self.assertIn('error_dict', info)