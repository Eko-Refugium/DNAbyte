from data_model import DNADS, RawData, EncodedData
from EncodingScheme_SimpleConstrained_mod import encoding_simpleConstrained

# create a DNADS instance with absolute paths
dnads_ins = DNADS(file_paths = ['/Users/fabian/Desktop/MI-DNA-DISC/tests/testfiles/Beethoven.jpg', 
                                '/Users/fabian/Desktop/MI-DNA-DISC/tests/testfiles/Ode_an_die_freude.odt'])

# alternative constructor
dnads_ins = DNADS.from_folder('../tests/testfiles/')
dnads_ins.info()    # call the info method

# create a RawData instance

#dna_ins = RawData.generate_bitstream(dnads_ins)
dnads_ins = RawData(dnads_ins)
dnads_ins.info()    # call the info method

# create an EncodedData instance

#dna_ins = EncodedData.generate_codewords(dnads_ins, simple_constrained)
dnads_ins = simple_constrained.encode(dnads_ins)

#EncodedData(dnads_ins, encoding_scheme = simple_constrained)
#dnads_ins.info()    # call the info method


# create a SynthesizedData instance
#dnads_ins = SynthesizedData(dnads_ins, synthesis_scheme = simple_constrained)

encoded_images_instance = EncodedData(dnads_ins, encoding_simpleConstrained, p=2, q=8, message_length=233)
