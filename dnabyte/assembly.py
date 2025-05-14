from scipy.constants import Avogadro
from itertools import groupby
import random
import numpy as np 

from dnabyte.oligo import Oligo, complement, translate_nested_list, translate_nested_list_poly,back_translate_nested_list_poly_binom_real, back_translate_nested_list_chain, back_translate_nested_list_poly, back_translate_nested_list, back_translate_nested_list_poly_binom,back_translate_nested_list_real,back_translate_nested_list_chain_real
from dnabyte.oligopool import OligoPool
from dnabyte.auxiliary import flatten_at_layer
from dnabyte.data import EncodedData, AssembledData
from dnabyte.ErrorChannels.synthesis_simulation import synthesis_simulation
import numpy as np


class SimulateAssembly:
    """
    Simulate the assembly of single stranded oligos into long chaines of double stranded DNA.
    """
    def __init__(self, params, logger=None):
            
            self.assembly_structure = params.assembly_structure
            self.encoding_scheme = params.encoding_scheme
            self.mean = params.mean
            self.std_dev = params.std_dev
            self.vol = params.vol
            self.hybridisation_steps = params.hybridisation_steps
            self.synthesis_method = params.synthesis_method
            self.theory = params.theory
            self.library = params.library

            # To allow simulation with large quantites of oligos, we will scale down the number, perform
            # the simulation, and then scale back up the results.
            scale = 1
            while self.mean > 1000:
                self.mean = self.mean / 10
                self.std_dev = self.std_dev / 10
                scale = scale * 10

            self.scale = scale

    def print_list_structure(self, lst, level=0):
        """
        Recursively prints the structure of a list of lists.

        Args:
            lst (list): The list to analyze.
            level (int): The current nesting level (used for indentation).
        """
        if isinstance(lst, list):
            print("  " * level + f"Level {level}: List with {len(lst)} elements")
            for item in lst:
                self.print_list_structure(item, level + 1)
        else:
            if isinstance(lst, Oligo):
                print("  " * level + f"Level {level}: {type(lst).__name__} ({lst})")
            elif isinstance(lst, OligoPool):
                for oligo in lst.pool:
                    print("  " * level + f"Level {level}: {type(oligo).__name__} ({oligo})")
                


    def process_tuple_list(tuple_list):
        def flip_tuple(t):
            return t[::-1]

        # Ensure the tuple containing 'A' is the first tuple
        if 'A' == tuple_list[1][0] or 'A' == tuple_list[1][-1]:
            tuple_list = (tuple_list[1], tuple_list[0])

        # Ensure 'A' is the first entry of the first tuple
        if tuple_list[0][0] != 'A':
            tuple_list = (flip_tuple(tuple_list[0]), flip_tuple(tuple_list[1]))

        return tuple_list
    

    # def add_empty_to_innermost(self, lst, library):
    #     if isinstance(lst, list):
    #         if all( isinstance(i, Oligo) for i in lst):
    #             testvar = [[Oligo([complement(lst[0].motifs[0],library.dictmotives), 'empty']), lst[0]], lst[1]]
    #             # print(testvar)
    #             # breakpoint()
    #             return [[Oligo([complement(lst[0].motifs[0],library.dictmotives), 'empty']), lst[0]], lst[1]]
    #         else:
    #             return [self.add_empty_to_innermost(i,library) for i in lst]
    #     return lst
    def add_empty_to_innermost(self, lst, library):
        finishedlist = []
        finishedlist.append(Oligo([complement(lst[0].motifs[0],library.dictmotives), 'empty']))
        finishedlist.extend(lst)
        return finishedlist

    def print_nested_list(self,nested_list, level=0):
        for item in nested_list:
            if isinstance(item, list):
                self.print_nested_list(item, level + 1)
            else:
                print(item,  'oligo')

    def nest_list(self, lst):
        """
        Recursively nests a list into the structure [[[[a, b], c], d], e].

        Args:
            lst (list): The input list to be nested.

        Returns:
            list: The nested list structure.
        """
        if len(lst) <= 1:
            return lst  # Base case: if the list has 1 or fewer elements, return it as is
        return [self.nest_list(lst[:-1]), lst[-1]]  # Recursively nest all but the last element with the last element


    def simulate(self, encoded_data):
    
        assembled_data = AssembledData(encoded_data)
        assembled_data.scale = self.scale

        if isinstance(encoded_data, EncodedData):

            if self.assembly_structure == 'synthesis':

                if self.synthesis_method == None:
                    assembled_data.data = encoded_data.data
                elif self.synthesis_method == 'nosynthpoly':
                    sequences_multiples = [seq for seq in encoded_data.data for _ in range(max(1, int(np.random.normal(self.mean, self.std_dev))))]
                    assembled_data.data = sequences_multiples
                else:
                    synthesized_data = synthesis_simulation(encoded_data.data, method_id=self.synthesis_method, mean=self.mean, std_dev=self.std_dev)
                    assembled_data.data = synthesized_data

            # linear assembly / linear encoding
            elif self.assembly_structure == 'linear_assembly' and  self.encoding_scheme == 'linear_encoding':
                
                restructuredlib = translate_nested_list(encoded_data.data, self.library.translationlibleft, self.library.translationlibright)
                
                if self.theory == 'yes':
                    poolingfinish = self.call_function_repeatedly_therory(restructuredlib, self.library, 'linear_encoding')
                    retranslated = back_translate_nested_list_chain(poolingfinish, self.library.translationlibleft, self.library.translationlibright)
                    assembled_data.data = retranslated

                if self.theory == 'no':
                    empyrestructuredlib = []
                    for i in range(len(restructuredlib)):
                        empyrestructuredlib.append(self.nest_list(self.add_empty_to_innermost(restructuredlib[i], self.library)))
                        
                    poolingfinish = self.call_function_repeatedly(empyrestructuredlib, self.library, 'linear_encoding')
                    poolingfinish = flatten_at_layer(poolingfinish, 1)
                    poolingfinish = OligoPool([]).join(pools=poolingfinish, mean=1, std_dev=0)
                    retranslated = back_translate_nested_list_chain_real(poolingfinish, self.library.translationlibleft, self.library.translationlibright)
                    assembled_data.data = retranslated

            # linear assembly / binomial encoding
            elif self.assembly_structure == 'linear_assembly' and  self.encoding_scheme == 'binomial_encoding':

                restructuredlib = translate_nested_list(encoded_data.data, self.library.translationlibleft, self.library.translationlibright)
                
                if self.theory == 'yes':
                    poolingfinish = self.call_function_repeatedly_therory(restructuredlib, self.library, 'binomial_encoding')
                    retranslated = back_translate_nested_list(poolingfinish, self.library.translationlibleft, self.library.translationlibright)
                    retranslated = flatten_at_layer(retranslated, 1)
                    assembled_data.data = retranslated

                elif self.theory == 'no':
                    poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'binomial_encoding')
                    retranslated = back_translate_nested_list_real(poolingfinish, self.library.translationlibleft, self.library.translationlibright)
                    assembled_data.data = retranslated
                    
            # positional assembly / linear encoding
            elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'linear_encoding':
                restructuredlib = translate_nested_list_poly(encoded_data.data, self.library)
                
                if self.theory == 'yes':
                    poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'linear_encoding')
                    retranslated = back_translate_nested_list_poly(poolingfinish, self.library)
                    assembled_data.data = retranslated
                if self.theory == 'no':
                    poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'linear_encoding')
                    retranslated = back_translate_nested_list_poly_binom_real(poolingfinish, self.library)
                    assembled_data.data = retranslated

            # positional assembly / binomial encoding
            elif self.assembly_structure == 'positional_assembly' and self.encoding_scheme == 'binomial_encoding':
                restructuredlib = translate_nested_list_poly(encoded_data.data, self.library)
                if self.theory == 'yes':
                    poolingfinish = self.call_function_repeatedly_therory(restructuredlib, self.library, 'binomial_encoding')
                    
                    retranslated = back_translate_nested_list_poly_binom(poolingfinish, self.library)
                    retranslated = flatten_at_layer(retranslated, 1)
                    assembled_data.data = retranslated
                elif self.theory == 'no':
                    poolingfinish = self.call_function_repeatedly(restructuredlib, self.library, 'binomial_encoding')
                    
                    retranslated = back_translate_nested_list_poly_binom_real(poolingfinish, self.library)
                    assembled_data.data = retranslated
            
            info = {}

            return assembled_data, info

        else: 
            raise ValueError("Data object is not supported.")


    def substitute_tuples_with_oligos(self, nested_list):
        for i in range(len(nested_list)):
            if isinstance(nested_list[i], tuple):
                nested_list[i] = Oligo(motifs=nested_list[i])
            elif isinstance(nested_list[i], list):
                nested_list[i] = self.substitute_tuples_with_oligos(nested_list[i])
            else:
                raise ValueError("Invalid data structure.")
        return nested_list


    def recursive_assembly(self, nested_list, mean, std_dev, hybridisation_steps):
        """
        Recursively assemble the data using the pool of oligos.
        """
        for item in nested_list:

            if all(isinstance(item, Oligo) for item in nested_list):
                combined_pool = OligoPool(oligo_list=nested_list, mean=self.mean, std_dev=self.std_dev)
                combined_pool.hybridise(hybridisation_steps)
                return combined_pool
        
            elif all(isinstance(item, OligoPool) or isinstance(item, Oligo) for item in nested_list):
                combined_pool = OligoPool.join(nested_list, mean=self.mean, std_dev=self.std_dev)
                combined_pool.hybridise(hybridisation_steps)
                return combined_pool
            
            else:
                for i in range(len(nested_list)):
                    if isinstance(nested_list[i], list):
                        nested_list[i] = self.recursive_assembly(nested_list[i], mean, std_dev, hybridisation_steps)

                combined_pool = OligoPool([]).join(pools=nested_list, mean=mean, std_dev=std_dev)
                combined_pool.hybridise(hybridisation_steps)
                return combined_pool
                

    def find_depth(self, nested_list):
        if not isinstance(nested_list, list):
            return 0
        elif not nested_list:
            return 1
        else:
            return 1 + max(self.find_depth(item) for item in nested_list)

    
    def call_function_repeatedly(self,nested_list,library,encoding):
        if encoding == 'binomial_encoding' and library.structure == 'positional_assembly':
            
            polypooling = self.recursive_assembly_theoretical_withmean(nested_list,library, mean=self.mean)
            
     
            times = self.find_depth(polypooling)
            for i in range(times):
                polypooling = self.recursive_assembly_messy_randompairing(polypooling,library)
               
           
        elif encoding == 'linear_encoding' and library.structure == 'positional_assembly':
            polypooling = self.recursive_assembly_theoretical_withmean(nested_list,library, mean=self.mean)
            
            times = self.find_depth(polypooling)

            for i in range(times):
                polypooling = self.recursive_assembly_messy_randompairing(polypooling,library)
           
        elif encoding == 'binomial_encoding' and library.structure == 'linear_assembly':
            times = self.find_depth(nested_list)

            polypooling = nested_list
            for i in range(times):
                polypooling = self.recursive_assembly_messy_randompairing(polypooling,library)
                
        elif encoding == 'linear_encoding' and library.structure == 'linear_assembly':
            
            times = self.find_depth(nested_list)
            polypooling = nested_list
            for i in range(times-1):
                polypooling = self.recursive_assembly_messy_randompairing(polypooling,library)
                
        return polypooling
    

    def recursive_assembly_messy_randompairing(self,nested_list,library):
        counter = 0
        
        if all(isinstance(subitem, Oligo) for subitem in nested_list) and len(nested_list)>1:
            combined_pool = OligoPool([]).join(pools=nested_list, mean=self.mean, std_dev=self.std_dev)
            
            tester=0
            checkup = []
            newoligos = []
            looping = True
            while looping:
                leftends = []
                rightends = []

                if newoligos != []:
                    combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)
                if  len(combined_pool.pool) == 0:
                    break

                newoligos = []
                for i in range(len(combined_pool.pool)):
                        
                    allinpool = combined_pool.pool[i]
                    allends = [allinpool.end('r','5'),allinpool.end('f','5'),allinpool.end('r','3'),allinpool.end('f','3')]
                    allvelidends = [item for item in allends if item is not None]

                    for item in allvelidends:
                        leftends.append([i,item])
                    newoligos.append(allinpool)
                
                leftends = sorted(leftends, key=lambda x: x[1])
                leftends = [list(group) for key, group in groupby(leftends, key=lambda x: x[1])]

                tester =+ 1
                motivrcounter = 0
                
                if len(checkup) == 0:
                    findingmotive = leftends[motivrcounter][0][1]
                    checkup.append(leftends[motivrcounter][0][1])
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1][:-1])
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]
                    else:
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                else:
                    while findingmotive in checkup:

                        findingmotive = leftends[motivrcounter][0][1]
                        motivrcounter+=1
                        
                        if motivrcounter == len(leftends)-1:
                            
                            looping = False
                            motivrcounter =-1
                            break
                    motivrcounter-=1
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1])  
                        checkup.append(leftends[motivrcounter][0][1][:-1]) 
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]  
                    else:
                        checkup.append(leftends[motivrcounter][0][1])
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                sublist1 = None
                sublist2 = None
                
                for sublist in leftends:
                    if all(item[1] == findingmotive for item in sublist):
                        sublist1 = sublist
                    elif all(item[1] == findingmotivecomp for item in sublist):
                        sublist2 = sublist   

                if sublist1 == None:
                    lensublist1=0
                else:
                    lensublist1 = len(sublist1)
                if sublist2 == None:
                    lensublist2=0
                else:
                    lensublist2 = len(sublist2)
                    random.shuffle(sublist2)
 
                newoligosupdated = []
                deletableoligos = []
                for j in range(min(lensublist1,lensublist2)):
                    
                    temp_pool = OligoPool([newoligos[sublist1[j][0]],newoligos[sublist2[j][0]]])
                    deletableoligos.append(newoligos[sublist1[j][0]])
                    deletableoligos.append(newoligos[sublist2[j][0]])
                    newoligotoadd = temp_pool.hybridise(1,library)
                    newoligos.append(newoligotoadd)
                
                elements_count = {element: deletableoligos.count(element) for element in set(deletableoligos)}
                
                result_list = []
                for item in newoligos:
                    if item in elements_count and elements_count[item] > 0:
                        elements_count[item] -= 1
                    else:
                        result_list.append(item)
                newoligos = result_list

            combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)
            return combined_pool
    
        elif all(isinstance(subitem, OligoPool) or isinstance(subitem, Oligo)  for subitem in nested_list)  and len(nested_list)>1:
            combined_pool = OligoPool([]).join(pools=nested_list, mean=self.mean, std_dev=self.std_dev)
            tester = 0
            checkup = []
            newoligos = []
            looping = True

            while looping:
                leftends = []
                rightends = []
                if newoligos != []:
                    combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)

                newoligos = []
                for i in range(len(combined_pool.pool)):
                    
                    if isinstance(combined_pool.pool[i], OligoPool):
                        for allinpool in combined_pool.pool[i].pool:
                            allends = [allinpool.end('r','5'),allinpool.end('f','5'),allinpool.end('r','3'),allinpool.end('f','3')]
                            allvelidends = [item for item in allends if item is not None]
                            
                            for item in allvelidends:
                                leftends.append([i,item])
                            
                            newoligos.append(allinpool)
                            
                        
                    else:
                        allinpool = combined_pool.pool[i]
                        allends = [allinpool.end('r','5'),allinpool.end('f','5'),allinpool.end('r','3'),allinpool.end('f','3')]
                        allvelidends = [item for item in allends if item is not None]
                        for item in allvelidends:
                            leftends.append([i,item])
                        newoligos.append(allinpool)
                            
                leftends = sorted(leftends, key=lambda x: x[1])
                leftends = [list(group) for key, group in groupby(leftends, key=lambda x: x[1])]
                
                tester =+ 1
                motivrcounter = 0
                
                if len(checkup) == 0:
                    findingmotive = leftends[motivrcounter][0][1]
                    checkup.append(leftends[motivrcounter][0][1])
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1][:-1])
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]
                    else:
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                else:
                    while findingmotive in checkup:

                        findingmotive = leftends[motivrcounter][0][1]
                        motivrcounter+=1
                        
                        if motivrcounter == len(leftends)-1:
                            
                            looping = False
                            motivrcounter =-1
                            break
                    motivrcounter-=1
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1])  
                        checkup.append(leftends[motivrcounter][0][1][:-1]) 
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]  
                    else:
                        checkup.append(leftends[motivrcounter][0][1])
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                sublist1 = None
                sublist2 = None
                
                for sublist in leftends:
                    if all(item[1] == findingmotive for item in sublist):
                        sublist1 = sublist
                    elif all(item[1] == findingmotivecomp for item in sublist):
                        sublist2 = sublist   

                if sublist1 == None:
                    lensublist1=0
                else:
                    lensublist1 = len(sublist1)
                if sublist2 == None:
                    lensublist2=0
                else:
                    lensublist2 = len(sublist2)
                    random.shuffle(sublist2)

                newoligosupdated = []
                deletableoligos = []
                for j in range(min(lensublist1,lensublist2)):
                    
                    temp_pool = OligoPool([newoligos[sublist1[j][0]],newoligos[sublist2[j][0]]])
                    deletableoligos.append(newoligos[sublist1[j][0]])
                    deletableoligos.append(newoligos[sublist2[j][0]])
                    newoligotoadd = temp_pool.hybridise(1,library)
                    newoligos.append(newoligotoadd)
                
                elements_count = {element: deletableoligos.count(element) for element in set(deletableoligos)}
                
                result_list = []
                for item in newoligos:
                    if item in elements_count and elements_count[item] > 0:
                        elements_count[item] -= 1
                    else:
                        result_list.append(item)
                newoligos = result_list
                
            combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)
            return combined_pool
        
        elif len(nested_list)==1:
            return nested_list[0]
        else:
            return [self.recursive_assembly_messy_randompairing(subelement,library) if isinstance(subelement, list) else subelement for subelement in nested_list]
    

    def recursive_assembly_messy(self,nested_list,library):
        counter=0
        
        if all(isinstance(subitem, Oligo) for subitem in nested_list) and len(nested_list)>1:
            combined_pool = OligoPool(oligo_list=nested_list, mean=self.mean, std_dev=self.std_dev)
            combined_pool.hybridise(self.hybridisation_steps,library)
            return combined_pool
    
        elif all(isinstance(subitem, OligoPool) or isinstance(subitem, Oligo)  for subitem in nested_list)  and len(nested_list)>1:
            combined_pool = OligoPool([]).join(pools=nested_list, mean=self.mean, std_dev=self.std_dev)
            combined_pool.hybridise(self.hybridisation_steps,library)
            return combined_pool
        
        elif len(nested_list)==1:
            return nested_list[0]
        else:

            return [self.recursive_assembly_messy(subelement,library) if isinstance(subelement, list) else subelement for subelement in nested_list]
        

    def call_function_repeatedly_therory(self,nested_list,library,encoding):
        """
        Calls the target function the specified number of times.

        :param target_function: The function to be called repeatedly.
        :param times: The number of times to call the function.
        :param args: Positional arguments to pass to the target function.
        :param kwargs: Keyword arguments to pass to the target function.
        """
        times=self.find_depth(nested_list)

        if encoding=='linear_encoding' and library.structure=="linear_assembly":
            for i in range(times-1):
            
                nested_list=self.recursive_assembly_chain_theoretical(nested_list,library)
            
        elif encoding=='binomial_encoding' and library.structure=="linear_assembly":
            for i in range(times-2):
            
                nested_list=self.recursive_assembly_theoretical(nested_list,library)

        elif (encoding=='linear_encoding' and library.structure=="positional_assembly"):
            for i in range(times-1):
            
                nested_list=self.recursive_assembly_theoretical(nested_list,library)

        elif (encoding=='binomial_encoding' and library.structure=="positional_assembly"):
            for i in range(times-2):
                nested_list=self.recursive_assembly_theoretical(nested_list,library)
       
        return nested_list


    def recursive_assembly_theoretical_withmean(self,nested_list,library,mean):
        counter=0

        if all(isinstance(subitem, Oligo) for subitem in nested_list) and len(nested_list)>1:
        
            alloligos = []
            for i in range(int(mean)):
                hybridized=nested_list[0]
                hybridizationpool=OligoPool([hybridized,nested_list[1]])
                hybridized=hybridizationpool.hybridise(1,library)

                for i in range(len(nested_list)-2):
                    hybridizationpool=OligoPool([]).join([hybridized,nested_list[i+2]],mean=1, std_dev=0)
                    
                    hybridized=hybridizationpool.hybridise(1,library)
                alloligos.append(hybridized)
            hybridizedpool = OligoPool(alloligos)

            return hybridizedpool
    
        elif all(isinstance(subitem, OligoPool) or isinstance(subitem, Oligo)  for subitem in nested_list)  and len(nested_list)>1:
            hybridized=nested_list[0]
        
            combined_pool = OligoPool([]).join([hybridized,nested_list[1]], mean=mean, std_dev=0)
            hybridized=combined_pool.hybridise(1,library)
            for i in range(len(nested_list)-2):
                hybridizationpool=OligoPool([]).join([hybridized,nested_list[i+2]], mean=mean, std_dev=0)
                hybridized=hybridizationpool.hybridise(1,library)
            return hybridized
        
        elif len(nested_list)==1:
            return nested_list[0]
        else:

            
            return [self.recursive_assembly_theoretical_withmean(subelement,library,mean) if isinstance(subelement, list) else subelement for subelement in nested_list]


    def recursive_assembly_theoretical(self,nested_list,library):
        counter=0
        
        if all(isinstance(subitem, Oligo) for subitem in nested_list) and len(nested_list)>1:
            
            hybridized=nested_list[0]
            hybridizationpool=OligoPool([hybridized,nested_list[1]])
            hybridized=hybridizationpool.hybridise(1,library)

            for i in range(len(nested_list)-2):
                hybridizationpool=OligoPool([]).join([hybridized,nested_list[i+2]], mean=1, std_dev=0)

                hybridized=hybridizationpool.hybridise(1,library)
                if counter==0:
                    
                    counter+=1
            counter2=0

            return hybridized
    
        elif all(isinstance(subitem, OligoPool) or isinstance(subitem, Oligo)  for subitem in nested_list)  and len(nested_list)>1:

            hybridized=nested_list[0]
        
            combined_pool = OligoPool([]).join([hybridized,nested_list[1]], mean=1, std_dev=0)
            hybridized=combined_pool.hybridise(1,library)
            for i in range(len(nested_list)-2):
                hybridizationpool=OligoPool([]).join([hybridized,nested_list[i+2]], mean=1, std_dev=0)
                hybridized=hybridizationpool.hybridise(1,library)
            return hybridized
        
        elif len(nested_list)==1:
            return nested_list[0]
        else:
            return [self.recursive_assembly_theoretical(subelement,library) if isinstance(subelement, list) else subelement for subelement in nested_list]
        

    def recursive_assembly_chain_theoretical(self,nested_list,library):
        counter=0
        
        if all(isinstance(subitem, Oligo) for subitem in nested_list):
            
            hybridized=nested_list[0]
            
            endoligo = Oligo([complement(nested_list[0].motifs[0],library.dictmotives),'empty'])

            hybridizedend=OligoPool([hybridized,endoligo])
            hybridized=hybridizedend.hybridise(1,library)
            hybridizationpool=OligoPool([]).join([hybridized,nested_list[1]], mean=1, std_dev=0)
            hybridized=hybridizationpool.hybridise(1,library)
            for i in range(len(nested_list)-2):
                hybridizationpool=OligoPool([]).join([hybridized,nested_list[i+2]], mean=1, std_dev=0)
                hybridized=hybridizationpool.hybridise(1,library)
                if counter==0:
                    
                    counter+=1
            
            return hybridized
    
        elif all(isinstance(subitem, OligoPool) or isinstance(subitem, Oligo)  for subitem in nested_list):
            hybridized=nested_list[0]
        
            combined_pool = OligoPool([]).join([hybridized,nested_list[1]], mean=1, std_dev=0)
            hybridized=combined_pool.hybridise(1,library)
            for i in range(len(nested_list)-2):
                hybridizationpool=OligoPool([]).join([hybridized,nested_list[i+2]], mean=1, std_dev=0)
                hybridized=hybridizationpool.hybridise(1,library)
            return hybridized
        else:
            return [self.recursive_assembly_chain_theoretical(subelement,library) if isinstance(subelement, list) else subelement for subelement in nested_list]
        

    def recursive_assembly_messy_randompairing_chain(self,nested_list,library):
        if all(isinstance(subitem, Oligo) for subitem in nested_list) and len(nested_list)>1:
            combined_pool = OligoPool([]).join(pools=nested_list, mean=self.mean, std_dev=self.std_dev)
            tester=0
            checkup = []
            newoligos = []
            looping = True
            while looping:
                leftends = []
                rightends = []
                if newoligos != []:
                    combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)

                newoligos = []
                for i in range(len(combined_pool.pool)):
                        
                    # else:
                    allinpool = combined_pool.pool[i]
                    allends = [allinpool.end('r','5'),allinpool.end('f','5'),allinpool.end('r','3'),allinpool.end('f','3')]
                    allvelidends = [item for item in allends if (item is not None)]
                    for item in allvelidends:
                        leftends.append([i,item])
                    newoligos.append(allinpool)
                                            
                leftends = sorted(leftends, key=lambda x: x[1])
                leftends = [list(group) for key, group in groupby(leftends, key=lambda x: x[1])]
                
                tester =+ 1
                motivrcounter = 0
                
                if len(checkup) == 0:
                    findingmotive = leftends[motivrcounter][0][1]
                    checkup.append(leftends[motivrcounter][0][1])
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1][:-1])
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]
                    else:
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                else:
                    while findingmotive in checkup:

                        findingmotive = leftends[motivrcounter][0][1]
                        motivrcounter+=1
                        
                        if motivrcounter == len(leftends)-1:
                            
                            looping = False
                            motivrcounter =-1
                            break
                    motivrcounter-=1
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1])  
                        checkup.append(leftends[motivrcounter][0][1][:-1]) 
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]  
                    else:
                        checkup.append(leftends[motivrcounter][0][1])
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                sublist1 = None
                sublist2 = None
                
                for sublist in leftends:
                    if all(item[1] == findingmotive for item in sublist):
                        sublist1 = sublist
                    elif all(item[1] == findingmotivecomp for item in sublist):
                        sublist2 = sublist   

                if sublist1 == None:
                    lensublist1=0
                else:
                    lensublist1 = len(sublist1)
                if sublist2 == None:
                    lensublist2=0
                else:
                    lensublist2 = len(sublist2)
                    random.shuffle(sublist2)

                newoligosupdated = []
                deletableoligos = []
                for j in range(min(lensublist1,lensublist2)):
                    
                    temp_pool = OligoPool([newoligos[sublist1[j][0]],newoligos[sublist2[j][0]]])
                    deletableoligos.append(newoligos[sublist1[j][0]])
                    deletableoligos.append(newoligos[sublist2[j][0]])
                    newoligotoadd = temp_pool.hybridise(1,library)
                    newoligos.append(newoligotoadd)
                
                elements_count = {element: deletableoligos.count(element) for element in set(deletableoligos)}
                
                result_list = []
                for item in newoligos:
                    if item in elements_count and elements_count[item] > 0:
                        elements_count[item] -= 1
                    else:
                        result_list.append(item)
                newoligos = result_list

            combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)
            return combined_pool
    
        elif all(isinstance(subitem, OligoPool) or isinstance(subitem, Oligo)  for subitem in nested_list)  and len(nested_list)>1:
            combined_pool = OligoPool([]).join(pools=nested_list, mean=self.mean, std_dev=self.std_dev)
            tester=0
            checkup = []
            newoligos = []
            looping = True
            while looping:
                leftends = []
                rightends = []
                if newoligos != []:
                    combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)

                newoligos = []
                for i in range(len(combined_pool.pool)):
                    
                    if isinstance(combined_pool.pool[i], OligoPool):
                        for allinpool in combined_pool.pool[i].pool:
                            allends = [allinpool.end('r','5'),allinpool.end('f','5'),allinpool.end('r','3'),allinpool.end('f','3')]
                            allvelidends = [item for item in allends if item is not None]
                            for item in allvelidends:
                                leftends.append([i,item])
                            
                            newoligos.append(allinpool)
                            
                        
                    else:
                        allinpool = combined_pool.pool[i]
                        allends = [allinpool.end('r','5'),allinpool.end('f','5'),allinpool.end('r','3'),allinpool.end('f','3')]

                        allvelidends = [item for item in allends if (item is not None)]
                        for item in allvelidends:
                            leftends.append([i,item])
                        newoligos.append(allinpool)
                        
                leftends = sorted(leftends, key=lambda x: x[1])
                leftends = [list(group) for key, group in groupby(leftends, key=lambda x: x[1])]
                tester =+ 1
                motivrcounter = 0
                
                if len(checkup) == 0:
                    findingmotive = leftends[motivrcounter][0][1]
                    checkup.append(leftends[motivrcounter][0][1])
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1][:-1])
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]
                    else:
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                else:
                    while findingmotive in checkup:

                        findingmotive = leftends[motivrcounter][0][1]
                        motivrcounter+=1
                        
                        if motivrcounter == len(leftends)-1:
                            
                            looping = False
                            motivrcounter =-1
                            break
                    motivrcounter-=1
                    if '*' in leftends[motivrcounter][0][1]:
                        checkup.append(leftends[motivrcounter][0][1])  
                        checkup.append(leftends[motivrcounter][0][1][:-1]) 
                        findingmotivecomp = leftends[motivrcounter][0][1][:-1]  
                    else:
                        checkup.append(leftends[motivrcounter][0][1])
                        checkup.append(leftends[motivrcounter][0][1]+'*')
                        findingmotivecomp = leftends[motivrcounter][0][1]+'*'
                    motivrcounter+=1
                sublist1 = None
                sublist2 = None
                
                for sublist in leftends:
                    if all(item[1] == findingmotive for item in sublist):
                        sublist1 = sublist
                    elif all(item[1] == findingmotivecomp for item in sublist):
                        sublist2 = sublist   

                if sublist1 == None:
                    lensublist1=0
                else:
                    lensublist1 = len(sublist1)
                if sublist2 == None:
                    lensublist2=0
                else:
                    lensublist2 = len(sublist2)
                    random.shuffle(sublist2)
                    
                newoligosupdated = []
                deletableoligos = []
                for j in range(min(lensublist1,lensublist2)):
                    
                    temp_pool = OligoPool([newoligos[sublist1[j][0]],newoligos[sublist2[j][0]]])
                    deletableoligos.append(newoligos[sublist1[j][0]])
                    deletableoligos.append(newoligos[sublist2[j][0]])
                    newoligotoadd = temp_pool.hybridise(1,library)
                    newoligos.append(newoligotoadd)
                
                elements_count = {element: deletableoligos.count(element) for element in set(deletableoligos)}
                
                result_list = []
                for item in newoligos:
                    if item in elements_count and elements_count[item] > 0:
                        elements_count[item] -= 1
                    else:
                        result_list.append(item)
                newoligos = result_list

            combined_pool = OligoPool([]).join(pools=newoligos, mean=1, std_dev=0)
            return combined_pool
        
        elif len(nested_list)==1:
            return nested_list[0]
        else:

            return [self.recursive_assembly_messy_randompairing_chain(subelement,library) if isinstance(subelement, list) else subelement for subelement in nested_list]
             







