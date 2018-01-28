"""
This is to make:

(1) a simple fasta file where a sequence is in a single line.
    Convert *.fasta extension to simple *.fa file

(2) and  prepend the accesion code in the header of each single line fasta
    For example: BWcode + '@'
    and save the resulting  .fa files

"""
import sys
import glob
from itertools import groupby

def fasta_iter(fh):
    """
    given a fasta file. yield tuples of header, sequence
    brentp @ https://www.biostars.org/p/710
    """
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next().strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

if __name__ == "__main__":
   
   for fa_path in glob.glob('*.fasta'): # Each fasta in current directoy
       BWcode = fa_path.split('.')[0].split('assembly_')[1]
       with open(fa_path, 'r') as inFasta, open(BWcode + '.fa', 'w') as outFas:
           for header, seq in fasta_iter(inFasta):
               tag = '>' + BWcode + '@' + header[1:]  
               outFas.write(tag + '\n')
               outFas.write(seq.strip() + '\n')
