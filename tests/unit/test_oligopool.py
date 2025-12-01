import unittest
from dnabyte.oligo import Oligo
from dnabyte.oligopool import OligoPool


class TestOligoPool(unittest.TestCase):
        
    def test_create_ds_oligo(self):
        """
        Test the creation of a single double stranded oligo by hybridising single stranded oligos in a pool
        """
        pool1 = OligoPool([Oligo(('a','A*')), Oligo(('a*','D'))])
        self.assertIsInstance(pool1, OligoPool)
        self.assertEqual(len(pool1.pool), 2)
        pool1.hybridise(1)
        self.assertEqual(len(pool1.pool), 1)
        
        self.assertIsInstance(pool1, OligoPool)
        self.assertIsInstance(pool1.pool[0], Oligo)
        self.assertIsInstance(str(pool1.pool[0]), str)


    def test_hybridise_double_stranded_oligos(self):
        """
        Test creation of a pool of double stranded oligos by simulating hybridisation of single stranded oligos
        that were pipettet into a single reaction chamber.
        """
        pool2 = OligoPool([Oligo(('a','A*')), Oligo(('a*','D'))], mean=100, std_dev=5)
        n_before = len(pool2.pool)
        pool2.hybridise(10000)
        n_after = len(pool2.pool)
        self.assertIsInstance(pool2, OligoPool)
        self.assertLessEqual(n_after, n_before)
        
    def test_hybridisation_of_three_oligos(self):
        """
        Test the creation and printing of the first tripple of oligos.
        """
        pool3 = OligoPool([Oligo(('b','B')), Oligo(('b*','A*')), Oligo(('a','B*'))], mean=100, std_dev=2)
        self.assertIsInstance(pool3, OligoPool)

        # How can one compare two OligoPools?
        #self.assertSetEqual(set(pool3.pool), set([Oligo(('b','B')), Oligo(('b*','A*')), Oligo(('a','B*'))]))

        self.assertEqual(len(set(pool3.pool)), 3)

        n_before = len(pool3.pool)
        pool3.hybridise(200000)
        n_after = len(pool3.pool)
        unique_pool = set(pool3.pool)

        self.assertGreaterEqual(len(unique_pool), 11)
        self.assertLessEqual(n_after, n_before)

    def test_hybridisation_of_nested_OligoPools(self):
        data_motifs = [
        [('b*','A*'), ('b','B'), ('a','B*'), # pool 1
        ('c*','A*'), ('c','C'), ('a','C*'), 
        ('d*','A*'), ('d','D'), ('a','D*')],

        [('a*','C*'), ('c','C'), ('c*','B*'), # pool 2
        ('a*','D*'), ('d','D'), ('d*','B*'), 
        ('a*','E*'), ('e','E'), ('e*','B*')],

        [('c', 'B'), ('c*', 'C'), ('b', 'C*'), # pool 3
        ('d', 'B'), ('d*', 'D'), ('b', 'D*'),
        ('e', 'B'), ('e*', 'E'), ('b', 'E*')]
        ]

        data_oligos = [[Oligo(motif) for motif in pool] for pool in data_motifs]
        data_pools = [OligoPool(pool, mean=100, std_dev=2) for pool in data_oligos]

        self.assertEqual(all(isinstance(pool, OligoPool) for pool in data_pools), True)
        data_pools_hybridised = [pool.hybridise(10000) for pool in data_pools]
        self.assertEqual(all(isinstance(pool, OligoPool) for pool in data_pools_hybridised), True)

        data_pools = OligoPool.from_oligo_pools(data_pools_hybridised)
        self.assertIsInstance(data_pools, OligoPool)

if __name__ == '__main__':
    unittest.main()




