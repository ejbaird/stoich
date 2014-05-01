#!/usr/bin/env python
#input: annotated reference genome and array.txt (tsv) of genes and differential expression between lines/ages
#output: chart of gene name, #C atoms, #N atoms, C/N ratio in both ancestor and evolved 

import argparse
import os 
import sys#, getopt
#import fileinput

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

def parse_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest='reference')
    parser.add_argument('-d', dest='diff_exp')
    args = parser.parse_args
    print args
    print args.reference, args.diff_exp
    return args

def AAelements():
    aa = {}
    aa[A] = {abbrev:'Ala', C:3, H:7, N:1, O:2}
    aa[R] = {abbrev:'Arg', C:6, H:14, N:4, O:2}
    aa[N] = {abbrev:'Asn', C:4, H:8, N:2, O:3}
    aa[D] = {abbrev:'Asp', C:7, H:7, N:1, O:4}
    aa[C] = {abbrev:'Cys', C:3, H:7, N:1, O:2, S:1}
    aa[Q] = {abbrev:'Gln', C:5, H:10, N:2, O:3}
    aa[E] = {abbrev:'Glu', C:5, H:9, N:1, O:4}
    aa[G] = {abbrev:'Gly', C:2, H:5, N:1, O:2}
    aa[H] = {abbrev:'His', C:6, H:9, N:3, O:2}
    aa[I] = {abbrev:'Ile', C:6, H:13, N:1, O:2}
    aa[L] = {abbrev:'Leu', C:6, H:13, N:1, O:2}
    aa[K] = {abbrev:'Lys', C:6, H:14, N:2, O:2}
    aa[M] = {abbrev:'Met', C:5, H:11, N:1, O:2, S:1}
    aa[F] = {abbrev:'Phe', C:9, H:11, N:1, O:2}
    aa[P] = {abbrev:'Pro', C:5, H:9, N:1, O:2}
    aa[S] = {abbrev:'Ser', C:3, H:7, N:1, O:3}
    aa[T] = {abbrev:'Thr', C:4, H:9, N:1, O:3}
    aa[W] = {abbrev:'Trp', C:11, H:12, N:2, O:2}
    aa[Y] = {abbrev:'Tyr', C:9, H:11, N:1, O:3}
    aa[V] = {abbrev:'Val', C:5, H:11, N:1, O:2}
    return aa

def main():
    args = parse_input()
    AA = AAelements() 

main()
