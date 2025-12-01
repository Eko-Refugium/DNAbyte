import unittest
from dnabyte.oligo import Oligo, complement, nucleotide_complement


class TestOligo(unittest.TestCase):
    def test_oligo(self):

        # test all ways to instantiate an Oligo object
        oligo1 = Oligo(('a', 'A*'))
        oligo2 = Oligo(sequence='ACTGTAGCTAGTCACTG')

        self.assertIsInstance(oligo1, Oligo)
        self.assertIsInstance(oligo2, Oligo)

        self.assertEqual(oligo1.type, 'single_stranded')
        self.assertEqual(oligo1.type, 'single_stranded')


    def test_str(self):

        # TEST A: create and print single stranded oligos 
        self.assertEqual(str(Oligo(('a','A*'))), "3':  (  a   -  A*  )  :5' \n")
        self.assertEqual(str(Oligo(sequence='ACTGTAGCTAGTCACTG')), "3': ACTGTAGCTAGTCACTG :5'")


    def test_complement(self):

        # TEST B: test nucleotide_complement
        self.assertEqual(nucleotide_complement('A'), 'T')
        self.assertEqual(nucleotide_complement('C'), 'G')
        self.assertEqual(nucleotide_complement('G'), 'C')
        self.assertEqual(nucleotide_complement('T'), 'A')
        self.assertEqual(nucleotide_complement(['A', 'C', 'G', 'T', 'A', 'A', 'A']), ['T', 'G', 'C', 'A', 'T', 'T', 'T'])

        # TEST C: test complement
        self.assertEqual(complement('a'), 'a*')
        self.assertEqual(complement(('b', 'C*')), ('b*', 'C'))




# Run the tests
if __name__ == '__main__':
    unittest.main()