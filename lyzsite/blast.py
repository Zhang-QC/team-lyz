from Bio.Blast.Applications import NcbiblastpCommandline 
from io import StringIO 
from Bio.Blast import NCBIXML 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord 
from Bio import SeqIO
import os

 # Create two sequence files
seq1 = SeqRecord(Seq("FQTWEEFSRAAEKLYLADPMKVRVVLKYRHVDGNLCIKVTDDLVCLVYRTDQAQDVKKIEKF"),
                   id="seq1")
seq2 = SeqRecord(Seq("FQTWEEFSRAEKLYLADPMKVRVVLRYRHVDGNLCIKVTDDLICLVYRTDQAQDVKKIEKF"),
                   id="seq2")
SeqIO.write(seq1, "seq1.fasta", "fasta")
SeqIO.write(seq2, "seq2.fasta", "fasta")


def write_seq(str1, str2, name1, name2):
	seq1 = SeqRecord(Seq("FQTWEEFSRAAEKLYLADPMKVRVVLKYRHVDGNLCIKVTDDLVCLVYRTDQAQDVKKIEKF"),
                   id="seq1")
	seq2 = SeqRecord(Seq("FQTWEEFSRAEKLYLADPMKVRVVLRYRHVDGNLCIKVTDDLICLVYRTDQAQDVKKIEKF"),
                   id="seq2")
	SeqIO.write(seq1, "seq1.fasta", "fasta")
SeqIO.write(seq2, "seq2.fasta", "fasta")


# Run BLAST and parse the output as XML
output = NcbiblastpCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5)()[0]
blast_result_record = NCBIXML.read(StringIO(output))

# Print some information on the result
for alignment in blast_result_record.alignments:
    for hsp in alignment.hsps:
        print('****Alignment****')
        print('sequence:', alignment.title)
        print('length:', alignment.length)
        print('e value:', hsp.expect)
        print(hsp.query)
        print(hsp.match)
        print(hsp.sbjct)

os.system(str(output))

