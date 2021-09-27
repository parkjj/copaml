# COPAML

This repository contains code accompanying the manuscript "Metaviromic identification of genetic hotspots of pathogenicity in SARS-CoV-2 using machine learning."

## Data acquisition
Coronavirus genome sequences were obtained from the Virus Pathogen Database and Analysis Resource (ViPR) database (https://www.viprbrc.org/brc/home.spg?decorator=vipr) and aligned with MAFFT version 7 (https://mafft.cbrc.jp/alignment/server/). A compressed version of the multiple sequence alignment (MSA) used in this study is uploaded here.

## Dependencies
The following software and module versions were used for training base machine learning models and obtaining performance scores:
```
Python  3.7.0
numpy 1.16.4   
pandas 0.24.2
scikit-learn 0.21.2  
biopython 1.76 
```
## Pipeline
There are three primary python scripts the perform the preprocessing, machine learning (ML) model training, and evaluation. Within each scripts, the working directory and input / output locations should be correctly defined.
* step1_parser.py
* step2_windowencoder.py
* step3_mltrain_[class type].py

To contain intermediate files generated throughout the process, the following empty directories should be initialized.
* windowseqs
* windowseqs_encoded
* windowseqs_scores_human
* windowseqs_scores_intersect
* windowseqs_scores_patho

### step1_parser.py
This script reads in the MSA and defines tiling 6bp windows with 1bp shifts across the alignment. It outputs the sample substrings for each window in the "windowseqs" directory.

### step2_windowencoder.py
This script reads in sample substrings for each window in the "windowseqs" directory and performs one-hot encoding, which is then exported to the  "windowseqs_encoded" directory. The script takes the window id as an argument, ex:
```
for x in {000000..100835}; do
    python step2_windowencoder.py $x
done
```

### step3_mltrain_[class type].py
This script reads in one-hot encodings for each window in the "windowseqs_encoded" directory, as well as the sample classification membership ("classprint_[class type].csv"), to train and evaluate the ML models. For this study, the following classifiers were used: Random Forests (RF), Support Vector Machines (SVM), Bernoulli Naïve Bayes (BNB), Gradient Boosting Classifiers (GBC) and Multi-layer Perceptron Classifiers (MLPC). The script takes the window id as an argument, ex:
```
for x in {000000..100835}; do
    python step3_mltrain_intersect.py $x
done
```

Class type ids are the following:
* "human" - corresponds to classification strategy A, which comprise coronavirus samples infecting human host.
* "patho" - corresponds to classification strategy B, which comprise all SARS-CoV-2, SARS-CoV, and MERS-CoV sample.
* "intersect" - corresponds to classification strategy C, which comprise SARS-CoV-2, SARS-CoV, and MERS-CoV samples infecting human hosts.


### Statistical meta-model
Scripts for calculating statistics for obtaining nucleotide-resolution coronavirus pathogenicity (COPA) scores are held in R scripts in the "stats" directory. Scores from the ML cross-valdiations performed above are pooled into a matrix where each row corresponds to a window and each column corresponds to a ML model and fold permutation. P-values, Q-values, and W-statistics are obtained by Wilcoxon rank sum test.


### Principal components analysis
Scripts for dimensionality reduction of encoded whole coronavirus genomes are held in R scripts in the "pca_analysis" directory. "pca_pre.R" converts the MSA to cell-based representations in a CSV file, "pca_run.R" performs both one hot encoding and PCA, "pca_post.R" creates visualizations with metadata labelling of the results from the PCA analysis.


## Contact
Please see manuscript for further details. For any questions, contact jonathan [dot] park [at] yale [dot] edu.

This study was conducted by the Sidi Chen lab at Yale University (https://sidichenlab.org/).
