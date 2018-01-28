# KASPTree

## About

KASPTree is a pipeline to make a phylogenetic tree using KASP marker sequences extracted from RenSeq (Resistance gene enrichment sequencing) data. The RenSeq protocol requires a bait library and that these KASP marker sequences be part of the bait library. 

Further details on this method can be found in the manuscript [http://biorxiv.org/cgi/content/short/248146v1](http://biorxiv.org/cgi/content/short/248146v1)

## Pre-requisites
### Python 2.7.13
We use Python programming ([https://www.python.org](https://www.python.org)) language version 2.7.13  from Anaconda ([https://anaconda.org](https://anaconda.org)) built (64-bit) with
Biopython [(http://biopython.org/)](http://biopython.org/)) library v1.68.

### Sequence quality trimming
We use Trimmomatic ([http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)) for read preprocessing and check the quality with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### *De novo* assembly software
We need assemblies of RenSeq data to fetch the KASP marker sequences for each accession. We suggest using CLC assembly cell ([https://www.qiagenbioinformatics.com/products/clc-assembly-cell/](https://www.qiagenbioinformatics.com/products/clc-assembly-cell/)) but other available programs for generating assemblies can also be used.

### BLAST+
For getting the best hit for each marker sequence, we used the BLAST+ command-line tools for alignment version 2.2.28 ([ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.28/))

## Steps

#### Step 1. Preprocess reads with Trimmomatic

#### Step 2. Generate assemblies with CLC assembly cell

#### Step 3. Generate single line fasta files (*.fa) and tag the contigs name with accession number as prefix 

Using python script: Single_Line_Fasta.py

#### Step 4. Make blast database for each assembly

```
accession1  makeblastdb -in accession1.fav  -out accession1 -dbtype nucl
accession2  makeblastdb -in accession2.fav  -out accession2 -dbtype nucl 
...
accessionN  makeblastdb -in accessionN.fav  -out accessionN -dbtype nucl
```

#### Step 5. BLAST KASP marker sequences to the assembly 
```
accession1  blastn -db accession1 -query markers.fasta  -out accession1.nofmt  -num_alignments 1  -num_descriptions 1 
...
accessionN  blastn -db accessionN -query markers.fasta  -out accessionN.nofmt  -num_alignments 1  -num_descriptions 1
 
```

#### Step 6. Parsing the blast output result using Biopython's NCBIStandalone and BlastParser 

Using Python script: `Parse_BLAST_output.py` 

#### Step 7. Generating the genotyping matrix and filter for markers that were missing in more than 60% of accessions and those with minor allele frequency less than 5 percent 

Using Python script: `Genotype_matrix.py`


#### Step 8. Construct UPGMA tree using Biopython (UPGMA-make-from-Trasposed.py script)

Using Python script: `UPGMA_tree.py`

#### Step 9. Visualization of tree on iTOL ([https://itol.embl.de/](https://itol.embl.de/))





 
