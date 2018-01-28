"""
 make UPGMA tree 

"""
import sys
import random
random.seed(123)
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.Consensus import *

#--------------- Transpose -------------------------------
def transPosed(csv_file = 'Table.tsv', transposedFasta = 'Table-T.fas'):
    #helper function
    def getVect(csv_file):
       d={}
       with open(csv_file, 'r') as inp:
         for line in inp:
           A=line.strip().split('\t')
           if len(A) >  2:
              d[A[0]]=A[1:]
       return d   

    #---------------------------------------------
    d=getVect(csv_file)
    ACC=d['marker']
    serd=[ s for s in d.keys() if s not in ['marker']]

    with open(transposedFasta,'w') as fas:
     for i, acc in enumerate(ACC):
       seq=''   
       for s in serd: 
         seq+=d[s][i] ## i is position
       fas.write('>'+ acc + '\n')
       fas.write(     seq + '\n')


#-------------- patching section for iTOL visualisation------------------------------
## Following are the parts from Biopython Consensus code 
from Bio.Phylo.NewickIO import Writer

def _get_comment(clade):
    if hasattr(clade, 'comment') and clade.comment:
        return _format_comment(str(clade.comment))
    else:
        return ''

def _info_factoryJ(self, plain, confidence_as_branch_length,
                      branch_length_only, max_confidence, format_confidence,
                      format_branch_length):
        """Return a function that creates a nicely formatted node tag."""
        if plain:
            # Plain tree only. That's easy.
            def make_info_string(clade, terminal=False):
                return _get_comment(clade)

        elif confidence_as_branch_length:
            # Support as branchlengths (eg. PAUP), ignore actual branchlengths
            def make_info_string(clade, terminal=False):
                if terminal:
                    # terminal branches have 100% support
                    return (':' + format_confidence % max_confidence) + _get_comment(clade)
                else:
                    return (':' + format_confidence % clade.confidence) + _get_comment(clade)

        elif branch_length_only:
            # write only branchlengths, ignore support
            def make_info_string(clade, terminal=False):
                return (':' + format_branch_length % clade.branch_length) + _get_comment(clade)

        else:
            # write support and branchlengths (e.g. .con tree of mrbayes)
            def make_info_string(clade, terminal=False):
                if (terminal or
                        not hasattr(clade, 'confidence') or
                        clade.confidence is None):
                    return (':' + format_branch_length
                            ) % (clade.branch_length or 0.0) + _get_comment(clade)
                else:
                    return (':' + format_branch_length + '[' + format_confidence + ']'   ## this  modify
                            ) % (clade.branch_length or 0.0, clade.confidence ) + _get_comment(clade) ## and this modify

        return make_info_string

#---------------------------------------
def bootstrap(msa, times=10):
    """Generate bootstrap replicates from a multiple sequence alignment object

    :Parameters:
        msa : MultipleSeqAlignment
            multiple sequence alignment to generate replicates.
        times : int
            number of bootstrap times.
    """

    length = len(msa[0]) 
    i = 0
    while i < times:
        i += 1
        item = None
        for j in range(length):
            col = random.randint(0, length - 1)
            if not item:
                item = msa[:, col:col + 1]
            else:
                item += msa[:, col:col + 1]
        yield item

def bootstrap_trees(msa, times, tree_constructor=None):
    """Generate bootstrap replicate trees from a multiple sequence alignment.

    :Parameters:
        msa : MultipleSeqAlignment
            multiple sequence alignment to generate replicates.
        times : int
            number of bootstrap times.
        tree_constructor : TreeConstructor
            tree constructor to be used to build trees.
    """

    msas = bootstrap(msa, times)
    for aln in msas:
        tree = tree_constructor.build_tree(aln)
        yield tree
        
## Until here from Biopython Consensus code
#-------------------------------------------------

if __name__ == "__main__":
   
  ## input Table.tsv
  transPosed(csv_file = 'Table.tsv', transposedFasta = 'Table-T.fas')

  aln = AlignIO.read(open('Table-T.fas'), 'fasta')
  dtc = DistanceTreeConstructor(DistanceCalculator("identity"), "upgma") #


  orig_tree = dtc.build_tree(aln) ## original Tree

  Phylo.write(orig_tree, "UPGMA-identity-Patched-orig.nwk", "newick") #

  bootrees = bootstrap_trees(aln, times = 100, tree_constructor = dtc)

  trees = list(bootrees)

  support_tree = get_support(orig_tree, trees)

  Writer._info_factory = _info_factoryJ

  Phylo.write(support_tree, "boot1K-UPGMA-identity.xml", "phyloxml")
  Phylo.write(support_tree, "boot1K-UPGMA-identity.nwk", "newick")

