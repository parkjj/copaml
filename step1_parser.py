import numpy as np
import sklearn as sklearn
from numpy import array
from numpy import argmax
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import re
import gc
import os
from Bio import SeqIO, AlignIO
from itertools import compress
import csv
import pandas as pd


os.chdir("your/directory/here")
fasta_sequences = AlignIO.read(open("data/msa.pir"),'fasta')

# Obtain kmers from sequences
def kmers(seq, size):
    return [seq[i:i+size].lower()
            for i in range(len(seq)-size+1)]

seqs = []
for i in range(0,len(fasta_sequences)):
    seqs.append(str(fasta_sequences[i].seq))

seqskmer = []
for i in range(0, len(seqs)):
    fullseq = seqs[i]
    fullseqkmer = kmers(fullseq, size=6)
    kmersub = fullseqkmer
    seqskmer.append(kmersub)

seqskmers = [list(x) for x in zip(*seqskmer)]

for i in range(0, len(seqskmers)):
   with open("windowseqs/window_" + str(i).zfill(6) + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(seqskmers[i])