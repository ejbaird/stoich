#!/usr/bin/env python
#input: annotated reference genome and array.txt (tsv) of genes and differential expression between lines/ages
#output: chart of gene name, #C atoms, #N atoms, C/N ratio in both ancestor and evolved 

import argparse
import os 
import sys
import fileinput

import config

from Bio import SeqIO
from Bio.seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from bumpy import random
import numpy as np
import operator
import random
import pandas as pd

def readinput():
    for line in file.input():
        process(line)
    return stuff?

def AAelements():
    aa = {}
    aa.A = {3let='Ala', C=3, H=7, N=1, O=2}
    aa.R = {3let='Arg', C=6, H=14, N=4, O=2}
    aa.N = {3let='Asn', C=4, H=8, N=2, O=3}
    aa.D = {3let='Asp', C=7, H=7, N=1, O=4}
    aa.C = {3let='Cys', C=3, H=7, N=1, O=2, S=1}
    aa.Q = {3let='Gln', C=5, H=10, N=2, O=3}
    aa.E = {3let='Glu', C=5, H=9, N=1, O=4}
    aa.G = {3let='Gly', C=2, H=5, N=1, O=2}
    aa.H = {3let='His', C=6, H=9, N=3, O=2}
    aa.I = {3let='Ile', C=6, H=13, N=1, O=2}
    aa.L = {3let='Leu', C=6, H=13, N=1, O=2}
    aa.K = {3let='Lys', C=6, H=14, N=2, O=2}
    aa.M = {3let='Met', C=5, H=11, N=1, O=2, S=1}
    aa.F = {3let='Phe', C=9, H=11, N=1, O=2}
    aa.P = {3let='Pro', C=5, H=9, N=1, O=2}
    aa.S = {3let='Ser', C=3, H=7, N=1, O=3}
    aa.T = {3let='Thr', C=4, H=9, N=1, O=3}
    aa.W = {3let='Trp', C=11, H=12, N=2, O=2}
    aa.Y = {3let='Tyr', C=9, H=11, N=1, O=3}
    aa.V = {3let='Val', C=5, H=11, N=1, O=2}
    return aa





