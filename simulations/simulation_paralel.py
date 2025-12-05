import os
import logging
import traceback
import time
from datetime import datetime
from tqdm import tqdm
from scipy.constants import Avogadro
import random
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed

from dnabyte.data_classes.base import Data
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.data_classes.insilicodna import InSilicoDNA

from dnabyte.binarization import Binarize
from dnabyte.encode import Encode
from dnabyte.library import Library
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.params import Params
from simulations.simulation_parameters import simulation_parameters


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

    def run(self, paralel=False, max_workers=4):
        # TODO: Implement parallel processing and allow the user to pass on arguments for parallel processing in the run method
        results = {}
        counter = 1

        if paralel:
            # Use ThreadPoolExecutor for parallel execution
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # Submit tasks to the executor
                future_to_params = {
                    executor.submit(self._run_single_simulation, params, counter): params
                    for counter, params in enumerate(self.parameters, start=1)
                }

                # Process completed tasks
                for future in as_completed(future_to_params):
                    params = future_to_params[future]
                    try:
                        result = future.result()
                        results.update(result)
                    except Exception as e:
                        self.simlogger.error('Error during parallel simulation for parameter set %s: %s', params.name, str(e))
                        self.simlogger.error(traceback.format_exc())

        else:
            # Sequential execution
            with tqdm(total=len(self.parameters), desc="Simulation Progress", unit="param") as pbar:
                for params in self.parameters:
                    name = params.name + '_' + str(counter)
                    counter += 1
                    results[name] = {}

                    if params.seed:
                        random.seed(params.seed + counter)

                    try:
                        self.simlogger.info('##################################################################################')
                        self.simlogger.info('SIMULATION SETTING: %s', params.name)
                        self.simlogger.info(params.__str__())

                        # Run the simulation for a single parameter set
                        single_result = self._run_single_simulation(params, counter)
                        results.update(single_result)

                    except Exception as e:
                        self.simlogger.error('Error during simulation for parameter set %s: %s', params.name, str(e))
                        self.simlogger.error(traceback.format_exc())
                        results[name]['status'] = 'FAILURE'
                        continue

                    # Update the progress bar
                    pbar.update(1)

        # Close the file handler
        self.handler.close()
        self.simlogger.removeHandler(self.handler)

        return results

    def _run_single_simulation(self, params, counter):
        """Helper function to run a single simulation."""
        results = {}
        name = params.name + '_' + str(counter)
        results[name] = {}

        if params.seed:
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
                data_obj = Data([params.filename])
                bin = Binarize(params)
                binary_code = bin.binarize(data_obj)


            except Exception as e:
                self.simlogger.info('STATUS: ERROR')
                self.simlogger.error('TYPE: %s', str(e))
                self.simlogger.error(traceback.format_exc())
                results[name]['status'] = 'FAILURE'
                return results

            # Generate the info
            duration = time.time() - start_time
            info = {
                'duration': duration,
                'length_of_bitsream': len(binary_code.data)
            }

            # Logging
            self.simlogger.info('STATUS: SUCCESS')
            self.simlogger.info('DURATION: %.2f seconds', duration)
            self.simlogger.info('LENGTH OF DATA: %d', len(binary_code.data))
            self.simlogger.info(binary_code.__str__())

            # Save to results
            results[name]['step1'] = info

#######################################################################################################################
##### STEP 2: ENCODE DATA #############################################################################################
#######################################################################################################################

            self.simlogger.info('STEP02: ENCODE DATA')
            start_time = time.time()

            try:
                #lib = Library(structure=params.assembly_structure, filename='./tests/testlibraries/' + params.library_name)
                enc = Encode(params, logger=self.simlogger)
                data_enc, info = enc.encode(binary_code)

            except Exception as e:
                self.simlogger.info('STATUS: ERROR')
                self.simlogger.error('TYPE: %s', str(e))
                self.simlogger.error(traceback.format_exc())
                results[name]['status'] = 'FAILURE'
                #continue

            # generate the info
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
                #continue

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
                #continue

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
                #continue

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

            self.simlogger.info('STEP06: PROCESS DATA')
            start_time = time.time()
            try:
                data_cor, info = enc.process(data_seq)

            except Exception as e:
                self.simlogger.info('STATUS: ERROR')
                self.simlogger.error('TYPE: %s', str(e))
                self.simlogger.error(traceback.format_exc())
                results[name]['status'] = 'FAILURE'
                #continue

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

            try:
                data_dec, valid, info = enc.decode(data_cor)

                if not valid:
                    self.simlogger.info('STATUS: ERROR')
                    self.simlogger.info('TYPE: decoded data does not match the original data')
                    raise ValueError("Decoding failed")

            except Exception as e:
                self.simlogger.info('STATUS: ERROR')
                self.simlogger.error('TYPE: %s', str(e))
                self.simlogger.error(traceback.format_exc())
                results[name]['status'] = 'FAILURE'
                #continue

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

            try:
                comparison, res = data_dec.compare(data_dec, binary_code, logger=self.simlogger)

            except Exception as e:
                self.simlogger.info('STATUS: ERROR')
                self.simlogger.error('TYPE: %s', str(e))
                self.simlogger.error(traceback.format_exc())
                results[name]['status'] = 'FAILURE'
                #continue

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
                bin.debinarize(data_dec, output_file_path='./tests/testfiles/testdecode')

            # try:
            #     RestoredData(data_dec, output_folder='./tests/testfiles/testdecode', job_identifier=self.job_identifier)

            except Exception as e:
                self.simlogger.info('STATUS: ERROR')
                self.simlogger.error('TYPE: %s', str(e))
                self.simlogger.error('TRACEBACK:' + traceback.format_exc())
                results[name]['status'] = 'FAILURE'
                #continue

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

        return results