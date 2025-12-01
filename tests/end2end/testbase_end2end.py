import os
import unittest
import logging
import traceback
import time
from datetime import datetime

from dnabyte.data_classes.base import Data
from dnabyte.data_classes.nucleobasecode import NucleobaseCode
from dnabyte.data_classes.binarycode import BinaryCode
from dnabyte.data_classes.insilicodna import InSilicoDNA

from dnabyte.binarize import Binarize
from dnabyte.encode import Encode
from dnabyte.library import Library
from dnabyte.synthesize import SimulateSynthesis
from dnabyte.store import SimulateStorage
from dnabyte.sequence import SimulateSequencing
from dnabyte.params import Params

class TestBase(unittest.TestCase):

    def __init__(self, methodName='runTest', params=None):
        super().__init__(methodName)
        self.params = params
        self.testlogger = None

    def setUp(self):
        # Set up logging configuration
        self.testlogger = logging.getLogger(__name__)
        self.testlogger.setLevel(logging.DEBUG)

    def tearDown(self):
        if hasattr(self, 'handler') and self.handler:
            self.handler.close()
            self.testlogger.removeHandler(self.handler)

    def test_logic(self):

        job_identifier = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_filename = os.path.join('tests', 'testlogs', f'job_{job_identifier}.log')
        self.handler = logging.FileHandler(log_filename)
        self.formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        self.handler.setFormatter(self.formatter)
        self.testlogger.addHandler(self.handler)

        try:
            self.testlogger.info('TEST SETTING: %s', self.params.name)
            self.testlogger.info(self.params.__str__())

#######################################################################################################################
##### STEP 1: BINARIZE DATA ###########################################################################################
#######################################################################################################################

            self.testlogger.info('STEP01: BINARIZE DATA')
            start_time = time.time()

            try:
                print('before data obj creation')
                data_obj = Data(['./tests/testfiles/' + self.params.filename])

                bin = Binarize(self.params)
                binary_code = bin.binarize(data_obj)

                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.info('DURATION: %.2f seconds', time.time() - start_time)
                self.testlogger.info('LENGTH OF DATA: %d', len(binary_code.data))
                self.testlogger.info(binary_code.__str__())
            except Exception as e:
                self.testlogger.info('STATUS: ERROR')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error('TRACEBACK:' + traceback.format_exc())

#######################################################################################################################
##### STEP 2: ENCODE DATA #############################################################################################
#######################################################################################################################

            self.testlogger.info('STEP02: ENCODE DATA')
            start_time = time.time()

            try:
                enc = Encode(self.params, logger=self.testlogger)
                data_enc, info = enc.encode(binary_code)
                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.info('DURATION: %.2f seconds', time.time() - start_time)
                self.testlogger.info('NUMBER OF CODEWORDS: %d', len(data_enc.data))
                self.testlogger.info('BARCODE LENGTH: %d', info['barcode_length'])
                self.testlogger.info(data_enc.__str__())
                # TODO: Add the info to the log
                #self.testlogger.info(info)

            except Exception as e:
                self.testlogger.info('STATUS: ERROR')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error(traceback.format_exc())
                self.fail(f"Encoding failed: {str(e)}")

#######################################################################################################################
##### STEP 3: SIMULATE SYNTHESIS ######################################################################################
#######################################################################################################################

            self.testlogger.info('STEP03: SIMULATE SYNTHESIS')
            start_time = time.time()

            try:
                syn = SimulateSynthesis(self.params, logger=self.testlogger)
                data_syn, info = syn.simulate(data_enc)
                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.info('DURATION: %.2f seconds', time.time() - start_time)
                self.testlogger.info('NUMBER OF CODEWORDS: %d', len(data_syn.data))
                self.testlogger.info(data_syn.__str__())

                # TODO: Add the info to the log
                #self.testlogger.info(info)

            except Exception as e:
                self.testlogger.info('STATUS: ERROR')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error(traceback.format_exc())
                self.fail(f"Synthesis simulation failed: {str(e)}")

#######################################################################################################################
##### STEP 4: SIMULATE STORAGE ########################################################################################
#######################################################################################################################

            self.testlogger.info('STEP04: SIMULATE STORAGE')
            start_time = time.time()
            try:
                sto = SimulateStorage(self.params, logger=self.testlogger)
                data_sto, info = sto.simulate(data_syn)
                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.info('DURATION: %.2f seconds', time.time() - start_time)
                self.testlogger.info('NUMBER OF STRAND BREAKS: %d', info['number_of_strand_breaks'])
                self.testlogger.info('NUMBER OF CODEWORDS: %d', len(data_sto.data))
                self.testlogger.info(data_sto.__str__())
            except Exception as e:
                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error(traceback.format_exc())
                self.fail(f"Storage simulation failed: %str(e)")

#######################################################################################################################
##### STEP 5: SIMULATE SEQUENCING #####################################################################################
#######################################################################################################################

            self.testlogger.info('STEP05: SIMULATE SEQUENCING')
            start_time = time.time()
            try:
                seq = SimulateSequencing(self.params, logger=self.testlogger)
                data_seq, info = seq.simulate(data_sto)
                data_seq = InSilicoDNA(data_seq.data)
                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.info('DURATION: %.2f seconds', time.time() - start_time)
                self.testlogger.info('NUMBER OF INTRODUCED ERRORS: %d', len(info))
                self.testlogger.info(data_seq.__str__())

            except Exception as e:
                self.testlogger.info('STATUS: ERROR')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error('Error during sequencing simulation: %s', str(e))
                self.testlogger.error(traceback.format_exc())
                self.fail(f"Sequencing simulation failed: %s")

#######################################################################################################################
##### STEP 6: PROCESSING ##############################################################################################
#######################################################################################################################

            self.testlogger.info('STEP06: PROCESS DATA')
            start_time = time.time()
            try:
                data_cor, info = enc.process(data_seq)
                data_cor = NucleobaseCode(data_cor.data)
                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.info('DURATION: %.2f seconds', time.time() - start_time)
                self.testlogger.info(data_cor.__str__())

                # TODO: Add the info to the log
                #self.testlogger.info(info)

            except Exception as e:
                self.testlogger.info('STATUS: ERROR')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error(traceback.format_exc())
                self.fail(f"Data processing failed: %s")

#######################################################################################################################
##### STEP 7: DECODING THE DATA #######################################################################################
#######################################################################################################################

            self.testlogger.info('STEP07: DECODE DATA')
            start_time = time.time()

            try:
                data_dec, valid, info = enc.decode(data_cor)

                # TODO: Add the info to the log
                #self.testlogger.info(info)

                if valid:
                    self.testlogger.info('STATUS: SUCCESS')
                    self.testlogger.info('DURATION: %.2f seconds', time.time() - start_time)
                    self.testlogger.info(data_dec.__str__())
                else:
                    self.testlogger.info('STATUS: ERROR')
                    self.testlogger.info('TYPE: decoded data does not match the original data')
                    self.fail("Decoding failed")
                    raise ValueError("Decoding failed")

            except Exception as e:
                self.testlogger.info('STATUS: ERROR')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error(traceback.format_exc())
                self.fail(f"Decoding failed: %s")

#######################################################################################################################
##### STEP 8: COMPARE DATA ############################################################################################
#######################################################################################################################
     
            self.testlogger.info('STEP08: COMPARE DATA')
            start_time = time.time()

            try:
                #print(data_dec.data, 'decoded data')
                #print(binary_code.data, 'raw data')
                comparison, res = data_dec.compare(data_dec, binary_code, logger=self.testlogger)

                if comparison == 'SUCCESS':
                    self.testlogger.info('STATUS: SUCCESS')
                    self.testlogger.info('DURATION: %.2f seconds' + "\n", time.time() - start_time)

                if comparison == 'ERROR':
                    self.testlogger.info('STATUS: ERROR')
                    self.testlogger.error('TRACEBACK:' + traceback.format_exc())
                    self.testlogger.error('COMPARISON RESULT: %s', res)
                    self.fail("Decoded data does not match the original data")

            except Exception as e:
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error(traceback.format_exc())
                raise


#######################################################################################################################
##### STEP 9: RECREATE ORIGINAL DATA ##################################################################################
#######################################################################################################################

            self.testlogger.info('STEP09: RESTORE DATA')
            start_time = time.time()

            try:

                bin.debinarize(data_dec, output_file_path='./tests/testfiles/testdecode')
                # RestoredData(data_dec, output_folder='./tests/testfiles/testdecode', job_identifier=job_identifier)
                self.testlogger.info('STATUS: SUCCESS')
                self.testlogger.info('DURATION: %.2f seconds' + "\n", time.time() - start_time)

            except Exception as e:
                self.testlogger.info('STATUS: ERROR')
                self.testlogger.error('TYPE: %s', str(e))
                self.testlogger.error('TRACEBACK:' + traceback.format_exc())

            # Assertions
            self.assertTrue(valid, "Decoding failed")
            self.assertEqual(comparison, "SUCCESS", "2Decoded data does not match the original data")

        except Exception as e:
            self.testlogger.error('Error during test: %s', str(e))
            self.testlogger.error(traceback.format_exc())
            self.fail(f"Test failed: {str(e)}")