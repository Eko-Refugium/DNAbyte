import random
from .oligo import Oligo, complement
import numpy as np

class OligoPool:
    """
        A class to represent a pool of oligonucleotides.

        There are two constructors for the OligoPool class:
        1. The first constructor takes a list of motifs and generates a pool of oligos with a random copy number for each motif.
        The distribution of the copy number is set to be normal, rounded to the next integer. The user can set the mean and
        standard deviation of the distribution.

        2. The second constructor takes a list of OligoPool objects and generates a new pool by combining the oligos in the
        pools of the input OligoPool objects. The user can also set the number of hybridisation events to be performed on the
        new pool.
        """

    def __init__(self, oligo_list, mean=1, std_dev=0):
        self.pool = []
        for oligo in oligo_list:
            self.pool.extend(self.pipette(oligo, mean, std_dev))

    def pipette(self, oligo, mean, std_dev):

        step = []
        num_oligos = int(round(np.random.normal(mean, std_dev)))
        if num_oligos <= 0:
            num_oligos = 1
        for _ in range(num_oligos):
            step.append(oligo)
        return step

    def join(self, pools, mean, std_dev):

        # Create a new OligoPool instance
        combined_pool = OligoPool([])

        oligo_pools = [arg for arg in pools if isinstance(arg, OligoPool)]
        oligos = [arg for arg in pools if isinstance(arg, Oligo)]

        # Extend the pool with elements from other OligoPool instances
        for pool in oligo_pools:
            combined_pool.pool.extend(pool.pool)

        # Use the pipette function to append a random copy of the oligos to the pool
        for oligo in oligos:
            combined_pool.pool.extend(self.pipette(oligo, mean, std_dev))
        
        return combined_pool


    @classmethod
    def from_oligo_pools(cls, oligo_pools):
        new_pool_motifs = []
        for pool in oligo_pools:
            #new_pool_motifs.extend(pool.oligo_list)
            new_pool_motifs.extend(pool.pool)
        
        new_pool = cls(new_pool_motifs)
        new_pool.pool = [oligo for pool in oligo_pools for oligo in pool.pool]
        
        return new_pool


    def __str__(self):
        if self.pool is not None:
            oligo_counts = {}
            
            # count occurences of all oligos
            
            for oligo in self.pool:
                oligo_key = str(oligo)
                oligo_key_mirrored = str(oligo.mirror())

                if (oligo_key in oligo_counts):
                    oligo_counts[oligo_key] += 1
                elif (oligo_key_mirrored in oligo_counts):
                    oligo_counts[oligo_key_mirrored] += 1
                else:
                    oligo_counts[oligo_key] = 1

            # sort oligos by count, descending
            sorted_oligos = dict(sorted(oligo_counts.items(), key=lambda item: item[1], reverse=True))

            # generate output string
            output = ""
            for oligo_key, count in sorted_oligos.items():
                output += f"{oligo_key}: {count} \n \n"
            return output

        else:
            return "Oligo pool is empty."


    def hybridise(self, n, library, info=False):

        # set default value for n
        if n is None:
            n = 10*len(self.pool)

        for _ in range(n):
            if len(self.pool) > 1:  # Ensure there are at least two oligos to hybridize
                # Randomly pick two distinct oligos
                oligo_A, oligo_B = random.sample(self.pool, 2)

                # Pair the selected oligos
                hybridised_oligo = self.pair(oligo_A, oligo_B,library)

                # Check if the hybridisation was successful
                if hybridised_oligo.motifs:
                    # Remove the original oligos from the pool
                    self.pool.remove(oligo_A)
                    self.pool.remove(oligo_B)

                    # Append the new hybridised oligo to the pool
                    self.pool.append(hybridised_oligo)

                    #if info:  # Optionally print information about the hybridisation event


        return self


    def pair(self, oligo_A, oligo_B,library):
        self.motif_dict=library.dictmotives
        
        if oligo_A.type == 'single_stranded' and oligo_B.type == 'single_stranded':  

            if complement(oligo_A.end('3'),self.motif_dict) == oligo_B.end('3'):
                new_motifs = ((None, oligo_A.motifs[0], oligo_A.motifs[1]), (oligo_B.motifs[1], oligo_B.motifs[0], None))

            elif complement(oligo_A.end('5'),self.motif_dict) == oligo_B.end('5'):
                new_motifs = ((oligo_A.motifs[0], oligo_A.motifs[1], None), (None, oligo_B.motifs[1], oligo_B.motifs[0]))                
            
            else:
                return Oligo(None)
            
            new_oligo = Oligo(new_motifs)
            return new_oligo


        if oligo_A.type == 'double_stranded' and oligo_B.type == 'double_stranded':  
            if oligo_A.end('f', '5') and oligo_B.end('r', '5') and oligo_A.end('f', '5') == complement(oligo_B.end('r', '5'),self.motif_dict):

                new_motifs = (oligo_A.motifs[0] + oligo_B.motifs[0][1:],    # forward strand
                                oligo_A.motifs[1][0:-1] + oligo_B.motifs[1])  # reverse strand
            
            # case: dsA forward 5' sticky end binds to dsB forward 5' sticky end
            elif oligo_A.end('f', '5') and oligo_B.end('f', '5') and oligo_A.end('f', '5') == complement(oligo_B.end('f', '5'),self.motif_dict):

                new_motifs = (oligo_A.motifs[0] + oligo_B.motifs[1][-2::-1],        # forward strand
                                oligo_A.motifs[1][0:-1] + oligo_B.motifs[0][::-1])    # reverse strand

            # case: the sticky end is on the reverse strand on the right and on the forward strand on the left
            elif oligo_A.end('r', '3') and oligo_B.end('f', '3') and oligo_A.end('r', '3') == complement(oligo_B.end('f', '3'),self.motif_dict):

                new_motifs = (oligo_A.motifs[0][0:-1] + oligo_B.motifs[0],    # forward strand
                                oligo_A.motifs[1] + oligo_B.motifs[1][1:])        # reverse strand
            
            # case: the sticky ends are on the reverse strands on the right
            
            elif oligo_A.end('r', '3') and oligo_B.end('r', '3') and oligo_A.end('r', '3') == complement(oligo_B.end('r', '3'),self.motif_dict):

                new_motifs = (oligo_A.motifs[0][:-1] + oligo_B.motifs[1][::-1],    # forward strand
                                oligo_A.motifs[1] + oligo_B.motifs[0][:-1][::-1])        # reverse strand

            # case: the sticky end is on the forward strand on the left and on the reverse strand on the right
            elif oligo_A.end('f', '3') and oligo_B.end('r', '3') and oligo_A.end('f', '3') == complement(oligo_B.end('r', '3'),self.motif_dict):

                new_motifs = ( oligo_B.motifs[0][0:-1] + oligo_A.motifs[0],    # forward strand
                                oligo_B.motifs[1] + oligo_A.motifs[1][1:])  # reverse strand

            # case: the sticky ends are on the forward strand on the left
            elif oligo_A.end('f', '3') and oligo_B.end('f', '3') and oligo_A.end('f', '3') == complement(oligo_B.end('f', '3'),self.motif_dict):

                new_motifs = ( oligo_B.motifs[1][1:][::-1] + oligo_A.motifs[0],    # forward strand
                                oligo_B.motifs[0][::-1] + oligo_A.motifs[1][1:])  # reverse strand

            # case: the sticky end is on the reverse strand on the left and the forward strand on the right
            elif oligo_A.end('r', '5') and oligo_B.end('f', '5') and oligo_A.end('r', '5') == complement(oligo_B.end('f', '5'),self.motif_dict):

                new_motifs = (oligo_B.motifs[0] + oligo_A.motifs[0][1:],    # forward strand
                                oligo_B.motifs[1][:-1] + oligo_A.motifs[1])  # reverse strand
            
            # case: the sticky ends are on the reverse strand on the left
            elif oligo_A.end('r', '5') and oligo_B.end('r', '5') and oligo_A.end('r', '5') == complement(oligo_B.end('r', '5'),self.motif_dict):

                new_motifs = (oligo_B.motifs[1][::-1] + oligo_A.motifs[0][1:],    # forward strand
                                oligo_B.motifs[0][:0:-1] + oligo_A.motifs[1])  # reverse strand

            else: 
                return Oligo(None) 
            
            return Oligo(new_motifs)
        
        else: # case where one oligo is single stranded and the other is double stranded

            if oligo_A.type == 'double_stranded' and oligo_B.type == 'single_stranded': 
                ds = oligo_A
                ss = oligo_B
            else:
                ds = oligo_B
                ss = oligo_A

            # case: double strand forward 5' sticky end binds to ss 5' sticky end
            if ds.end('f', '5') and ds.end('f', '5') == complement(ss.end('5'),self.motif_dict):

                new_motifs=(ds.motifs[0] + (None,),                                 # forward strand
                                ds.motifs[1][:-1] + (ss.motifs[1], ss.motifs[0]))   # reverse strand
            
            # case: ds reverse 3' sticky end binds to ss 3' sticky end
            elif ds.end('r', '3') and ds.end('r', '3') == complement(ss.end('3'),self.motif_dict):

                new_motifs = (ds.motifs[0][0:-1] + (ss.motifs[0], ss.motifs[1]),    # forward strand
                                ds.motifs[1] + (None,))                             # reverse strand
            
            # case: ds forward 3' sticky end binds to ss 5' sticky end
            elif ds.end('f', '3') and ds.end('f', '3') == complement(ss.end('3'),self.motif_dict):  

                new_motifs = ((None,) + ds.motifs[0],                               # forward strand
                                (ss.motifs[1], ss.motifs[0]) + ds.motifs[1][1:])    # reverse strand
            
            # case: ds reverse 5' sticky end binds to ss 5' sticky end
            elif ds.end('r', '5') and ds.end('r', '5') == complement(ss.end('5'),self.motif_dict):

                new_motifs = ((ss.motifs[0], ss.motifs[1]) + ds.motifs[0][1:],      # forward strand
                                (None,) + ds.motifs[1])                             # reverse strand

            else:
                return Oligo(None)
            
            return Oligo(new_motifs)

