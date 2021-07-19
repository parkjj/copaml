import sys
import numpy as np
import sklearn as sklearn
from numpy import array
from numpy import argmax
import re
import gc
import os
from itertools import compress
import csv
import pandas as pd
import timeit
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import BernoulliNB
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold


os.chdir("your/directory/here")
input1 = sys.argv[1]
print(input1)
idclass = pd.read_csv('classprint_patho.csv', header=None)
idclass = idclass.transpose().values.tolist()
idclass = idclass[0]

windowseq = pd.read_csv('windowseqs_encoded/window_'+str(input1)+'_encoded.csv', index_col=0)

X = windowseq
y = idclass
X = array(X)
y = array(y)
state = 17

clf_svc = SVC(kernel='linear', C=1, random_state=state)
clf_rf = RandomForestClassifier(random_state=state, n_estimators=100)
clf_bnb = BernoulliNB()
clf_mlpc = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=state)
clf_gbc = GradientBoostingClassifier(random_state=state)

scores_svc = cross_val_score(clf_svc, X, y, cv=5)
scores_rf  = cross_val_score(clf_rf, X, y, cv=5)
scores_gnb = cross_val_score(clf_gnb, X, y, cv=5)
scores_mlpc = cross_val_score(clf_mlpc, X, y, cv=5)
scores_gbc  = cross_val_score(clf_gbc, X, y, cv=5)

scores_df = pd.DataFrame([scores_svc, scores_rf, scores_gnb,scores_mlpc,scores_gbc], index=['svc', 'rf', 'gnb', 'mlpc', 'gbc'])
scores_df.to_csv('windowseqs_scores_patho/window_'+str(input1)+'_scores.csv')
