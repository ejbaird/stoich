stoich
======

CTurner's project for determining counts and ratios of elements in amino acids in an annotated genome

Type the following into the command line:

    python elstoich.py -g <name of genome reference file> -d <name of differential expression file>

For example, 

    python elstoich.py -g REL606.2.gbk -d arrays.txt

###Output of program

gene name | #C atoms | #N atoms | C:N ratio   
(for all genes in the given file)

###Structure of program

1. [done] input differential expression array and reference
2. [done] determine atom counts for each amino acid
3. [done] create Python-based genbank dictionary (genename, locustag, translation)
4. create list of genes from diff-exp-array (important to do?)
5. determine #C, #N, and C:N per gene for the genome

###Reference paper

Cooper, T., Rozen, D., and R. Lenski. 2003. Parallel changes in gene expression after 20,000 generations of evolution in _E. coli._ PNAS. 100(3):1072-1077. 
