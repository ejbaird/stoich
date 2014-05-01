stoich
======

CTurner's project for determining counts and ratios of elements in amino acids in an annotated genome

Type the following into the command line:

    python elstoich.py -g <name of genome reference file> -d <name of differential expression file>

For example, 

    python elstoich.py -g REL606.2.gbk -d arrays.txt

###Structure of program

1. [done] input differential expression array and reference
2. [done] determine atom counts for each amino acid
3. create Python-based genbank dictionary, esp. with gene name and translation
4. create list of genes from diff_exp_array (important?)
5. determine #C, #N, and C:N per gene for the genome
