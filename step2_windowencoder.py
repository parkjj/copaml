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
import sys

os.chdir("your/directory/here")

input1 = sys.argv[1]
print(input1)

# String helper functions
def string_to_array(input_string):
    input_string = input_string.lower()
    input_string = re.sub('[\W]', 'z', input_string)
    input_string = re.sub('[^acgtnz]', 'n', input_string)
    output_array = np.array(list(input_string))
    return output_array

# Label encoder
label_enc = LabelEncoder()
label_enc.fit(np.array(['a','c','g','t','z','n']))

# One hot encoder
def one_hot(my_array):
    int_enc= label_enc.transform(my_array)
    oh_enc = OneHotEncoder(sparse=False, dtype=int, categories=[range(6)])
    int_enc = int_enc.reshape(len(int_enc), 1)
    oh_enc = oh_enc.fit_transform(int_enc)
    oh_enc = np.delete(oh_enc, -3, 1) # This deletes the 3rd column from the right - which encode "N"
    return oh_enc

windowseq = pd.read_csv('windowseqs/window_'+str(input1)+'.csv', header=None)
windowseq = windowseq.transpose().values.tolist()

seqsencoded = []
for i in range(0, len(windowseq)):
    seqsencoded.append(one_hot(string_to_array(windowseq[i][0])).ravel().tolist())

seqsencoded_df = pd.DataFrame(seqsencoded)
seqsencoded_df.to_csv('windowseqs_encoded/window_'+str(input1)+'_encoded.csv')