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

class AminoAcid(object):

    def __init__(self, letter, name, C, H, N, O, S=0):
        self.letter = letter
        self.name = name
        self.C = C
        self.H = H
        self.N = N
        self.O = O
        self.S = S

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
            mytuple = (gene_name, locus_tag, translation)
            geneinfo.append(mytuple)
    print '\n', "Genome gene list (1st 5): ", '\n', geneinfo[:5]
    print len(geneinfo)
    return geneinfo

def geneCN(sequence, AAdict):
    '''goal: take a single gene's AA sequence and obtain from it #C, #N, and C:N'''
    for AA in sequence:
        cn_tuple = ('#C', '#N', 'C:N')

def allgenesCN(genedict):
    CN_dict = {}
    for gene in genedict:
        try:
            run geneCN for each genedict[ #references a dict of tuples
        except KeyError:
            pass #?

def main():
    args = parse_input()
    AA = getAminoAcids() 
    rec = make_record(args.genref)
    geneAA = get_geneAAs(rec)
    CNcounts = allgeneCN(geneAA, AA)

main()
