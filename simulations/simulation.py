import os
import logging
import traceback
import time
from datetime import datetime
from tqdm import tqdm
from scipy.constants import Avogadro
import random
import pickle

# TODO: Check imports and remove unused ones
from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.data_classes.insilicodna import InSilicoDNA

from dnabyte.encode import Encode
from dnabyte.library import Library
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.params import Params
from dnabyte.binarize import Binarize

class Simulation():

    def __init__(self, simulation_parameters, debug=False):
        self.simlogger = logging.getLogger(__name__)
        self.simlogger.setLevel(logging.DEBUG)
        self.parameters = simulation_parameters

        self.job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_filename = os.path.join('simulations', 'simlogs', f'job_{self.job_identifier}.log')
        self.handler = logging.FileHandler(log_filename)
        self.formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        self.handler.setFormatter(self.formatter)
        self.simlogger.addHandler(self.handler)

    def run(self, paralel=False):

        # TODO: Implement parallel processing and allow the user to pass on arguments for paralel processing in the run method
        if paralel:
            pass

        results = {}
        counter = 1

        with tqdm(total=len(self.parameters), desc="Simulation Progress", unit="param") as pbar:
            for params in self.parameters:
                name = params.name + '_' + str(counter)
                counter += 1
                results[name] = {}

                if hasattr(params, 'seed') and params.seed:
                    random.seed(params.seed + counter)

                try:
                    self.simlogger.info('##################################################################################')
                    self.simlogger.info('SIMULATION SETTING: %s', params.name)
                    self.simlogger.info(params.__str__())

    #######################################################################################################################
    ##### STEP 1: BINARIZE DATA ###########################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP01: BINARIZE DATA')
                    start_time = time.time()

                    try:
                        bin = Binarize(params)

                        data_obj = Data(file_paths=params.file_paths)
                        binary_code = bin.binarize(data_obj)

                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    # generate the info
                    duration = time.time() - start_time
                    info = {
                        'duration': duration,
                        'length_of_bitsream': len(binary_code.data)
                    }

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds', duration)
                    self.simlogger.info('LENGTH OF DATA: %d', len(binary_code.data))
                    self.simlogger.info(binary_code.__str__())

                    # save to results
                    results[name]['step1'] = info

    #######################################################################################################################
    ##### STEP 2: ENCODE DATA #############################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP02: ENCODE DATA')
                    start_time = time.time()

                    try:
                        coder = Encode(params, logger=self.simlogger)
                        data_enc, info = coder.encode(binary_code)

                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    # generate the info
                    #print('info', info)
                    duration = time.time() - start_time
                    info['duration'] = duration
                    
                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds', duration)
                    self.simlogger.info('NUMBER OF CODEWORDS: %d', len(data_enc.data))
                    self.simlogger.info('BARCODE LENGTH: %d', info['barcode_length'])
                    self.simlogger.info(data_enc.__str__())

                    # save to results
                    results[name]['step2'] = info

    #######################################################################################################################
    ##### STEP 3: SIMULATE SYNTHESIS ######################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP03: SIMULATE SYNTHESIS')
                    start_time = time.time()

                    try:
                        syn = SimulateSynthesis(params, logger=self.simlogger)
                        data_syn, info = syn.simulate(data_enc)
                        
                        
                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue
                    # print('data_syn', data_syn.data)
                    # generate the info
                    duration = time.time() - start_time
                    info['duration'] = duration

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds', duration)
                    self.simlogger.info('NUMBER OF CODEWORDS: %d', len(data_syn.data))
                    self.simlogger.info(data_syn.__str__())

                    # save to results
                    results[name]['step3'] = info

    #######################################################################################################################
    ##### STEP 4: SIMULATE STORAGE ########################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP04: SIMULATE STORAGE')
                    start_time = time.time()
                    try:
                        sto = SimulateStorage(params, logger=self.simlogger)
                        data_sto, info = sto.simulate(data_syn)
                        
                    except Exception as e:
                        self.simlogger.info('STATUS: SUCCESS')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    # generate the info
                    duration = time.time() - start_time
                    info['duration'] = duration

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds', duration)
                    self.simlogger.info('NUMBER OF STRAND BREAKS: %d', info['number_of_strand_breaks'])
                    self.simlogger.info('NUMBER OF CODEWORDS: %d', len(data_sto.data))
                    self.simlogger.info(data_sto.__str__())

                    # save to results
                    results[name]['step4'] = info

    #######################################################################################################################
    ##### STEP 5: SIMULATE SEQUENCING #####################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP05: SIMULATE SEQUENCING')
                    start_time = time.time()
                    try:
                        seq = SimulateSequencing(params, logger=self.simlogger)
                        data_seq, info = seq.simulate(data_sto)
                        data_seq = InSilicoDNA(data_seq.data)

                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error('Error during sequencing simulation: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    # generate the info
                    duration = time.time() - start_time
                    info['duration'] = duration

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds', duration)
                    self.simlogger.info('NUMBER OF INTRODUCED ERRORS: %d', len(info))
                    self.simlogger.info(data_seq.__str__())

                    # save to results
                    results[name]['step5'] = info

    #######################################################################################################################
    ##### STEP 6: PROCESSING ##############################################################################################
    #######################################################################################################################
                    # print('data_seq', data_seq.data)
                    self.simlogger.info('STEP06: PROCESS DATA')
                    start_time = time.time()
                    try:
                        data_cor, info = coder.process(data_seq)
                    
                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue
                    # print('data_cor', data_cor.data)
                    # generate the info
                    duration = time.time() - start_time
                    info['duration'] = duration

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds', duration)
                    self.simlogger.info('NUMBER OF CODEWORDS: %d', len(data_cor.data))
                    self.simlogger.info(data_cor.__str__())

                    # save to results
                    results[name]['step6'] = info

    #######################################################################################################################
    ##### STEP 7: DECODING THE DATA #######################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP07: DECODE DATA')
                    start_time = time.time()
                    # print('data_cor', data_cor.data)
                    try:
                        data_dec, valid, info = coder.decode(data_cor)

                        if not valid:
                            self.simlogger.info('STATUS: ERROR')
                            self.simlogger.info('TYPE: decoded data does not match the original data')
                            raise ValueError("Decoding failed")

                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    # generate the info
                    duration = time.time() - start_time
                    info['duration'] = duration

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds', duration)
                    self.simlogger.info(data_dec.__str__())

                    # save to results
                    results[name]['step7'] = info

    #######################################################################################################################
    ##### STEP 8: COMPARE DATA ############################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP08: COMPARE DATA')
                    start_time = time.time()

                    # print('data_dec', data_dec.data)
                    # print('binary_code', binary_code.data)
                    # print('len data_dec', len(data_dec.data))
                    # print('len binary_code', len(binary_code.data))

                    try:
                        comparison, res = data_dec.compare(data_dec, binary_code, logger=self.simlogger)

                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    if comparison == 'ERROR':
                        self.simlogger.info('STATUS: ERROR')
                        raise ValueError("Decoded data does not match the original data")

                    # generate the info
                    duration = time.time() - start_time
                    info = {}
                    info['duration'] = duration
                    info['comparison'] = res

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds' + "\n", duration)

                    # save to results
                    results[name]['step8'] = res

    #######################################################################################################################
    ##### STEP 9: RECREATE ORIGINAL DATA ##################################################################################
    #######################################################################################################################

                    self.simlogger.info('STEP09: RESTORE DATA')
                    start_time = time.time()

                    try:
                        success = bin.debinarize(data=binary_code, output_directory='./simulations/simfiles/')

                        #RestoredData(data_dec, output_folder='./simulations/simfiles/', job_identifier=self.job_identifier)

                    except Exception as e:
                        self.simlogger.info('STATUS: ERROR')
                        self.simlogger.error('TYPE: %s', str(e))
                        self.simlogger.error('TRACEBACK:' + traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    # generate the info
                    duration = time.time() - start_time
                    info = {
                        'duration': duration
                    }

                    # logging
                    self.simlogger.info('STATUS: SUCCESS')
                    self.simlogger.info('DURATION: %.2f seconds' + "\n", duration)

                    # save to results
                    results[name]['step9'] = info

                    results[name]['status'] = 'SUCCESS'


                except Exception as e:
                    self.simlogger.error('Error during simulation: %s', str(e))
                    self.simlogger.error(traceback.format_exc())
                    results[name]['status'] = 'FAILURE'
                    continue

                # Update the progress bar
                finally: 
                    pbar.update(1)

            # Close the file handler
            self.handler.close()
            self.simlogger.removeHandler(self.handler)

            # Save the results to a pickle file
            pickle_filename = os.path.join('simulations', 'simlogs', f'res_{self.job_identifier}.pickle')
            try:
                with open(pickle_filename, 'wb') as pickle_file:
                    pickle.dump(results, pickle_file)
                self.simlogger.info(f"Results saved to {pickle_filename}")
            except Exception as e:
                self.simlogger.error(f"Failed to save results to {pickle_filename}: {str(e)}")
                self.simlogger.error(traceback.format_exc())

            return results