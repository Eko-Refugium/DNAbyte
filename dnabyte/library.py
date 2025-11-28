import csv

class Library:
    """
    A class to represent a DNA oligo library.
    """

    def __init__(self, structure='linear_assembly', library=None, **kwargs):
        """
        Initialize a Library object.

        :param structure: The structure of the library indicated by a string 'simple' (default) or 'poly'.
        :param library: The library data, default is None.
        :param kwargs: Additional keyword arguments.
        """
        self.structure = structure
        
        if library:
            self.library = library
        elif 'filename' in kwargs and self.structure == 'linear_assembly':
            self.library = self.read_library(kwargs['filename'])
            self.leftmotives, self.rightmotives, self.dictmotives, self.translationlibleft, self.translationlibright = self.motive_pairs()
        
            
        elif 'filename' in kwargs and self.structure == 'positional_assembly':
            self.messages, self.generic, self.position  = self.read_library_poly(kwargs['filename'])
            self.library = [self.messages, self.generic, self.position]
            self.messageleft, self.messageright, self.genericlib, self.connectorlib, self.dictmotives = self.motive_pairs_poly()


    @classmethod
    def read_library(cls, file_path):
        """
        Read the library data from a file.

        :param file_path: Path to the file containing the library data.
        :return: A list of DNA sequences.
        """
        # TODO: Fallunterscheidung mit fasta und .csv
        # TODO: make it an actual class method
        # TODO: check if read in library is correct

        
        DNAs = []
        with open(file_path) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:    
                DNAs.append(row[0])
        return DNAs      
    
    @classmethod
    def read_library_poly(cls, file_path):
        
        # TODO: check if read in library is correct
        
        DNAs = []
        with open(file_path) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:    
                DNAs.append(row[0])

        messages = DNAs[1:DNAs.index('Generic')]
        generic = DNAs[DNAs.index('Generic')+1:DNAs.index('Connector')]
        positions = DNAs[DNAs.index('Connector')+1:]

        return messages,  generic, positions
    
    def complementmap(self, string):
        compliment = ''
        for i in string:
            if i == 'A':
                compliment += 'T'
            if i == 'T':
                compliment += 'A'
            if i == 'C':
                compliment += 'G'
            if i == 'G':
                compliment += 'C'
        compliment = compliment[::-1]
        return compliment
    
    def motive_pairs(self):
        """
        Generate motive pairs from the library data.

        :param library: The library data.
        :return: A tuple containing left motifs, right motifs, a dictionary of motives,
             translation library for left motifs, and translation library for right motifs.
        """
        DNAs = self.library
        oligolen = len(DNAs[0])
        oligolenhalf = oligolen // 2
        leftmotifs = []
        rightmotifs = []

        for dna in DNAs:
            left_motif = dna[:oligolenhalf]
            right_motif = dna[oligolenhalf:]

            if left_motif not in leftmotifs and self.complementmap(left_motif) not in leftmotifs:
                leftmotifs.append(left_motif)
            if right_motif not in rightmotifs and self.complementmap(right_motif) not in rightmotifs:
                rightmotifs.append(right_motif) 

        leftmotifs.sort()
        rightmotifs.sort()

        translationlibleft = {}
        translationlibright = {}
        dictmotives = {}

        for i, motif in enumerate(leftmotifs):
            translationlibleft[motif] = chr(97 + i)
            translationlibleft[self.complementmap(motif)] = chr(97 + i) + '*'
            dictmotives[chr(97 + i)] = chr(97 + i) + '*'
            dictmotives[chr(97 + i) + '*'] = chr(97 + i)

        for i, motif in enumerate(rightmotifs):
            translationlibright[motif] = chr(65 + i)
            translationlibright[self.complementmap(motif)] = chr(65 + i) + '*'
            dictmotives[chr(65 + i)] = chr(65 + i) + '*'
            dictmotives[chr(65 + i) + '*'] = chr(65 + i)
        
        dictmotives['empty'] = 'empty*'
        dictmotives['empty*'] = 'empty'
        dictmotives[None] = None

        return leftmotifs, rightmotifs, dictmotives, translationlibleft, translationlibright
        
    def motive_pairs_poly(self):
        """
        Generate motive pairs from the library data.

        :param library: The library data.
        :return: A tuple containing left motifs, right motifs, a dictionary of motives,
             translation library for left motifs, and translation library for right motifs.
        """

        messagelibleft = {}
        message_libright = {}
        genericlib = {}
        connectorlib = {}
        dictmotives = {}
        messagesleft =[]
        messagesright = []

        for i in range(len(self.messages)):
            messagesleft.append(self.messages[i][:len(self.messages[i])//2])
            messagesright.append(self.messages[i][len(self.messages[i])//2:])
        # messagesleft = list(set(messagesleft))
        # messagesright = list(set(messagesright))
        # messagesleft.sort()
        # messagesright.sort()

        for i, motif in enumerate(messagesleft):
            messagelibleft[motif] = 'ml' + str(i)
            messagelibleft[self.complementmap(motif)] = 'ml' + str(i) + '*'
            dictmotives['ml' + str(i)] = 'ml' + str(i) + '*'
            dictmotives['ml' + str(i) + '*'] = 'ml' + str(i)

        for i, motif in enumerate(messagesright):
            message_libright[motif] = 'mr' + str(i)
            message_libright[self.complementmap(motif)] = 'mr' + str(i) + '*'
            dictmotives['mr' + str(i)] = 'mr' + str(i) + '*'
            dictmotives['mr' + str(i) + '*'] = 'mr' + str(i)

        for i, motif in enumerate(self.generic):
            genericlib[motif] = 'g' + str(i)
            genericlib[self.complementmap(motif)] = 'g' + str(i) + '*'
            dictmotives['g' + str(i)] = 'g' + str(i) + '*'
            dictmotives['g' + str(i) + '*'] = 'g' + str(i)

        for i, motif in enumerate(self.position):
            connectorlib[motif] = 'c' + str(i)
            connectorlib[self.complementmap(motif)] = 'c' + str(i)+ '*'
            dictmotives['c' + str(i)] = 'c' + str(i) + '*'
            dictmotives['c' + str(i) + '*'] = 'c' + str(i)

        dictmotives['empty'] = 'empty*'
        dictmotives['empty*'] = 'empty'
        dictmotives[None]= None

        return messagelibleft, message_libright, genericlib, connectorlib, dictmotives
        
