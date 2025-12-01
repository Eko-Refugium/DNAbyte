
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

sequences = ["AGGGGC",
             "AGGGC",
             "AGGGGGC",
             "AGGAGC",
             "AGGGGG"]

records = (SeqRecord(Seq(s)) for s in sequences)

SeqIO.write(records, "./tests/msa_example.fasta", "fasta")

# alignment = AlignIO.read("./tests/msa_example.fasta", "fasta")
# align_array = np.array([list(rec) for rec in alignment], np.character)



# for record in alignment:
#     print("%s - %s" % (record.seq, record.id))

# print(alignment)
# print(align_array)

cline = MuscleCommandline(input="./tests/msa_example.fasta", out="./tests/msa.txt")
print(cline)







