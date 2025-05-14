import tarfile
import os 
import gzip
import random


class DNADS:
    """
    DNADS is the base class for all data types in the MI-DNA-DISC project.

    It can be instantiated 
     - with a list of absolute paths to files or, alternatively
     - with a path to a folder using the alternative constructor from_folder.

    :param file_paths: A list of absolute paths to files.    
    """
    def __init__(self, file_paths):
        for file_path in file_paths:
            if not os.path.isfile(file_path):
                raise ValueError(f"File path {file_path} does not point to an existing file.")
        self.file_paths = file_paths
        #self.size = self.calculate_total_bytes()

    @classmethod
    def from_folder(cls, folder_path):
        if not os.path.isdir(folder_path):
            raise ValueError(f"Folder path {folder_path} does not point to an existing directory.")
        file_paths = [os.path.join(folder_path, filename) for filename in os.listdir(folder_path)]
        return cls(file_paths)

    def calculate_total_bytes(self):
        total_bytes = 0
        for file_path in self.file_paths:
            total_bytes += os.path.getsize(file_path)
        return total_bytes

    def __str__(self):
        output = "Type: " + str(type(self)) + "\n"
        #output += "Bytes:" + str(self.size) + "\n"
        output += "File paths: " + str(self.file_paths) + "\n"
        return output

    def binarize(self):
        """
        Converts the files in the DNADS object to a binary sequence.

        :return: A string of binary data.
        """
        # Create a compressed tar archive file
        folder_path = os.path.dirname(self.file_paths[0])
        archive_name = folder_path + '/archive.tar.gz'

        with tarfile.open(archive_name, 'w:gz') as archive:
            for file_path in self.file_paths:
                archive.add(file_path)
        
        size = os.path.getsize(archive_name)

        with open(archive_name, 'rb') as file:
            binary_data = file.read()
            binary_sequence = ''.join(format(byte, '08b') for byte in binary_data)

        os.remove(archive_name)

        obj = RawData(self)
        obj.size = size
        obj.file_paths = self.file_paths
        obj.data = binary_sequence
        obj.length = len(binary_sequence)

        return obj



class RawData(DNADS):
    """
    RawData is a subclass of DNADS.

    It consists of a compressed bitstream of the files in the DNADS object. For testing purposes, the RawData object can also directly be created from a bitstream.
    
    :param dnads: A DNADS object.
    """
    def init(self, bitstream):
        """
        Creates a RawData object from a bitstream.

        :param bitstream: A string of binary data.
        """
        self.file_paths = None
        instance = self.__new__(self)
        instance.data = bitstream
        instance.length = len(bitstream)
        instance.size = len(bitstream) // 8  # Size in bytes
        instance.file_paths = None

        return instance

    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "Type: " + str(type(self)) + "\n"
        # output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output


class EncodedData(DNADS):
    """
    EncodedData is a subclass of DNADS.

    Its main attribute is the data object, which consists of a multi-leveled list. At the highest level, data consists of a list of codewords, there is one list for every codeword in the data object. At the lowest level there are lists of motifs or nucleotides, that are to be hybridised in a single pool. The list structure represents the hierarchical structure of the pooling process. 

    :param raw_data: A RawData object.
    :param encoding_scheme: An EncodingScheme object.
    """
    #TODO: make it an actual subclass of DNADS 
    def __init__(self, raw_data, **kwargs):
        #self.file_paths = raw_data.file_paths
        #self.size = raw_data.size
        self.data = raw_data
        if "file_paths" in kwargs:
            self.file_paths = kwargs['file_paths']


    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "Type: " + str(type(self)) + "\n"
        # #output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "Number of Codewords : " + str(len(self.data)) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output


class AssembledData(DNADS):
    """
    AssembledData is a subclass of DNADS.

    It consists of a list of codewords generated from the list of codewords of the EncodedData object. This step simulates the process of synthesizing DNA from the codewords. In this process the codewords are generated in stochastic number of copies and errors are introduced.

    :param encoded_data: An EncodedData object.
    :param synthesis_simulator: A SynthesisSimulator object.
    """

    def __init__(self, encoded_data=None):
        if encoded_data is not None and isinstance(encoded_data, EncodedData):
            self.file_paths = encoded_data.file_paths
            #self.size = encoded_data.size

        self.data = None        

    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "\n" + "Type: " + str(type(self)) + "\n"
        # #output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output


class StoredData(DNADS):
    """
    StoredData is a subclass of DNADS. It consists of a list of stored codewords generated from the list of synthesized codewords of the SynthesizedData object.

    :param synthesized_data: A SynthesizedData object.
    :param storage_simulator: A StorageSimulator object.
    """
    def __init__(self, assembled_data=None):
        if assembled_data is not None and isinstance(assembled_data, AssembledData):                                        
            self.file_paths = assembled_data.file_paths
            self.size = assembled_data.size
            self.data = assembled_data.data
            self.scale = assembled_data.scale
        else:
            self.file_paths = None
            self.size = None
            self.data = assembled_data
            self.scale = 1


    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "\n" + "Type: " + str(type(self)) + "\n"
        # if self.file_paths is not None:
        #     output += "File paths: " + str(self.file_paths) + "\n"
        #     output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output


class SequencedData(DNADS):
    """
    SequencedData is a subclass of DNADS. It consists of a list of codewords generated from the list of stored codewords of the StoredData object.

    :param data: A StoredData object.
    :param sequencing_simulator: A SequencingSimulator object.
    """
    def __init__(self, data):
        if data is not None and isinstance(data, StoredData):
            self.file_paths = data.file_paths
            self.size = data.size
            self.data = data.data
            self.scale = data.scale
        else:
            self.file_paths = None
            self.size = None
            self.data = data
            self.scale = 1
        

    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "\n" + "Type: " + str(type(self)) + "\n"
        # if self.file_paths is not None:
        #     output += "File paths: " + str(self.file_paths) + "\n"
        #     output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output

class CorrectedData(DNADS):

    def __init__(self, sequenced_data):
        if sequenced_data is not None and isinstance(sequenced_data, SequencedData):
            self.file_paths = sequenced_data.file_paths
            self.size = sequenced_data.size
            self.data = sequenced_data.data
        else:
            self.file_paths = None
            self.size = None
            self.data = sequenced_data

    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "\n" + "Type: " + str(type(self)) + "\n"
        # if self.file_paths is not None:
        #     output += "File paths: " + str(self.file_paths) + "\n"
        #     output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output


class DecodedData(DNADS):
    """
    DecodedData is a subclass of DNADS. It consists of a bitstream created from the list of sequenced codewords of the SequencedData object.

    :param sequenced_data: A SequencedData object.
    :param encoding_scheme: An EncodingScheme object.
    """
    def __init__(self, sequenced_data):
        self.file_paths = sequenced_data.file_paths
        self.size = sequenced_data.size
        self.data = None
        
    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "\n" + "Type: " + str(type(self)) + "\n"
        # if self.file_paths is not None:
        #     output += "File paths: " + str(self.file_paths) + "\n"
        #     output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output

    
# class RestoredData(DNADS):
#     """
#     RestoredData is a subclass of DNADS.

#     It consists of a list of files created from the bitstream of the DecodedData object by decompressing the tar archive given in form of a binary sequence.

#     :param decoded_data: A DecodedData object.
#     """
#     def __init__(self, decoded_data, output_folder):
        
#         self.file_paths = decoded_data.file_paths
#         bitstream = decoded_data.data

#         binary_data = bytes([int(bitstream[i:i+8], 2) for i in range(0, len(bitstream), 8)])

#         # Get the directory of this script
#         script_dir = os.path.dirname(os.path.realpath(__file__))

#         # Make output_folder relative to the script's directory
#         output_folder = os.path.join(script_dir, output_folder)

#         # Ensure the output_folder directory exists
#         os.makedirs(output_folder, exist_ok=True)

#         # Save the binary data to a temporary tar archive file
#         with open(os.path.join(output_folder, 'temp_archive.tar.gz'), 'wb') as file:
#             file.write(binary_data)
    
#         # Extract files from the tar archive
#         with tarfile.open(os.path.join(output_folder, 'temp_archive.tar.gz'), 'r:gz') as archive:
#             for member in archive.getmembers():
#                 # Ignore the path in the archive and extract the file to output_folder
#                 member.name = os.path.basename(member.name)
#                 archive.extract(member, output_folder)

#         self.status = "SUCCESS"

class RestoredData:
    """
    RestoredData is a class that reconstructs the original files from a binary sequence.

    :param raw_data: A RawData object.
    :param output_folder: The folder where the restored files will be saved.
    """
    def __init__(self, raw_data, output_folder, job_identifier):
        """
        Reconstructs the original files from a binary sequence.

        :param raw_data: A RawData object.
        :param output_folder: The folder where the restored files will be saved.
        """
        # Suppress the specific DeprecationWarning
        # Warnings.filterwarnings("ignore", category=DeprecationWarning)

        self.file_paths = raw_data.file_paths
        self.data = raw_data.data

        # Convert binary sequence back to bytes
        byte_data = bytes(int(self.data[i:i+8], 2) for i in range(0, len(self.data), 8))

        # Create a temporary tar archive file
        archive_name = os.path.join(output_folder, job_identifier + '_restored.tar.gz')
        with open(archive_name, 'wb') as file:
            file.write(byte_data)

        # Extract the tar archive
        with tarfile.open(archive_name, 'r:gz') as archive:
            for member in archive.getmembers():
                # Remove the path information and add "_restored" to the file name
                original_name = os.path.basename(member.name)
                restored_name = f"{os.path.splitext(original_name)[0]}_restored{os.path.splitext(original_name)[1]}"
                member.name = restored_name
                archive.extract(member, output_folder)

        # Remove the temporary tar archive file
        #os.remove(archive_name)


    def __str__(self):
        output = "DATA: " + str(self.data)[:50] + "...\n"

        # output = "\n" + "Type: " + str(type(self)) + "\n"
        # if self.file_paths is not None:
        #     output += "File paths: " + str(self.file_paths) + "\n"
        #     output += "Bytes:" + str(self.calculate_total_bytes()) + "\n"
        # output += "DATA: " + str(self.data)[:200] + "...\n"
        return output
    
    



