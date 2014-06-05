#!/usr/bin/env python
#input: annotated reference genome and array.txt (tsv) of genes and differential expression between lines/ages
#output: chart of gene name, #C atoms, #N atoms, C/N ratio in both ancestor and evolved 

import argparse
import os 
import sys

from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', dest='genref')
    parser.add_argument('-d', dest='diff_exp')
    # parse_args() is a function; needs parens
    args = parser.parse_args()
    return args

class Gene(object):  '''creating new data structure. a class describes a way of storing data and methods and memory (how to store in memory) that are accessible from same entry point (i.e., Gene). class describes how to store data. when you instantiate an instance (i.e., create an object), that's one output of that class (class is blueprint, object is what it makes). Gene, here, inherits from object. everything is of type object (most basic type of object). self is the instancw that the method is called on'''
    def __init__(self, gene_name, locus_tag, sequence):
        self.gene_name = gene_name
        self.locus_tag = locus_tag
        self.sequence = sequence
        self.N = 0
        self.C = 0

class AminoAcid(object):
    
    def __init__(self, letter, name, C, H, N, O, S=0):
        self.letter = letter
        self.name = name
        self.C = C
        self.H = H
        self.N = N
        self.O = O
        self.S = S
    
    def __str__(self):
        return '{} ({}): C={} N={} H={} O={}'.format(self.letter, self.name, self.C, self.N, self.H, self.O)

def getAminoAcids():
    aa = {}
    aa['A'] = AminoAcid('A', 'Ala', 3, 7, 1, 2)
    aa['R'] = AminoAcid('R', 'Arg', 6, 14, 4, 2)
    aa['N'] = AminoAcid('N', 'Asn', 4, 8, 2, 3) 
    aa['D'] = AminoAcid('D', 'Asp', 7, 7, 1, 4)
    aa['C'] = AminoAcid('C', 'Cys', 3, 7, 1, 2, S=1)
    aa['Q'] = AminoAcid('Q', 'Gln', 5, 10, 2, 3)
    aa['E'] = AminoAcid('E', 'Glu', 4, 9, 1, 4)
    aa['G'] = AminoAcid('G', 'Gly', 2, 5, 1, 2)
    aa['H'] = AminoAcid('H', 'His', 6, 9, 3, 2)
    aa['I'] = AminoAcid('I', 'Ile', 6, 13, 1, 2)
    aa['L'] = AminoAcid('L', 'Leu', 6, 13, 1, 2)
    aa['K'] = AminoAcid('K', 'Lys', 6, 14, 2, 2)
    aa['M'] = AminoAcid('M', 'Met', 5, 11, 1, 2, S=1)
    aa['F'] = AminoAcid('F', 'Phe', 9, 11, 1, 2)
    aa['P'] = AminoAcid('P', 'Pro', 5, 9, 1, 2)
    aa['S'] = AminoAcid('S', 'Ser', 3, 7, 1, 3)
    aa['T'] = AminoAcid('T', 'Thr', 4, 9, 1, 3)
    aa['W'] = AminoAcid('W', 'Trp', 11, 12, 2, 2)
    aa['Y'] = AminoAcid('Y', 'Tyr', 9, 11, 1, 3)
    aa['V'] = AminoAcid('V', 'Val', 5, 11, 1, 2)
    return aa

def make_record(genome):
    '''obtain BioPython object of the refseq'''

    it = SeqIO.parse(genome, "genbank")
    record = it.next()
    return record

def get_geneAAs(record): 
    geneinfo = []
    for i in record.features:
        if i.type == 'CDS':
            locus_tag = i.qualifiers['locus_tag'][0]
            try:
                gene_name = i.qualifiers['gene'][0]
            except KeyError:
                gene_name = 'none'
            try:
                translation = i.qualifiers['translation'][0] 
            except KeyError:
                translation = 'none'
            geneinfo.append(Gene(gene_name, locus_tag, translation))
    print '\n', "Genome gene list (1st 5): ", '\n', geneinfo[:5]
    print len(geneinfo)
    return geneinfo

def geneCN(sequence, AAdict):
    '''goal: take a single gene's AA sequence and obtain from it #C, #N, and C:N'''
    numC = 0
    numN = 0
    for AA in sequence:
        numC += AAdict[AA].C
        numN += AAdict[AA].N
    return numC, numN

def allgenesCN(genelist, AAdict):
 #       try:
 #           run geneCN for each genelist[ #references a dict of tuples
 #       except KeyError:
 #           pass #?

    for gene in genelist:
        numC, numN = geneCN(gene.sequence, AAdict)
        gene.C = numC
        gene.N = numN

def main():
    args = parse_input()
    AA = getAminoAcids() 
    rec = make_record(args.genref)
    geneAA = get_geneAAs(rec)
    CNcounts = allgeneCN(geneAA, AA)

main()
