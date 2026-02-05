#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 16:22:25 2025

@author: jesusav
"""


'''
Lbraries
'''
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.metrics import ConfusionMatrixDisplay, RocCurveDisplay,PrecisionRecallDisplay
from sklearn.model_selection import cross_val_score
from sklearn import tree
from sklearn.metrics import accuracy_score, confusion_matrix


"""Lee. RF functions"""
sys.path.append('/RF_testing/Scripts/')
from DMI_RF_functions import preprocessing_dataset
from DMI_RF_functions import split_fit_rf
from DMI_RF_functions import make_confusion_matrix
from DMI_RF_functions import make_ROC_curve
from DMI_RF_functions import make_precision_recall_curve
from DMI_RF_functions import make_feature_importance_plot



'''
Input tables
'''
elm_classes = pd.read_csv("/RF_testing/Datasets/elm_classes_taxons_Jun2024.csv")
human_classes = elm_classes[elm_classes['ELM_Hs']==1]['ELMIdentifier']

PRS = pd.read_csv('/RF_testing/Datasets/PRS_20210413_featrs.tsv', sep='\t', index_col= 0)
RRS_1 = pd.read_csv('/RF_testing/Datasets/RRSv4_1_20210428_featrs.tsv', sep='\t', index_col= 0)


all_features = ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein1', 'DomainFreqinProteome1']
all_features_renamed = ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']
exclude_feature= []
# 


'''
Model rerun
'''
all_features= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein1', 'DomainFreqinProteome1']
all_features_renamed= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']
exclude_feature= []

df_1, X_1, y_1= preprocessing_dataset(PRS, RRS_1, all_features, all_features_renamed)

exclude_feature= []
outputs = []
for boot in range(10): #10 random splits
    output_1= split_fit_rf(X_1, y_1)
    outputs.append(output_1)

#make_confusion_matrix(outputs)
make_ROC_curve(outputs, 25)
make_precision_recall_curve(outputs, fontsize=25)
make_feature_importance_plot(outputs, all_features_renamed, fontsize=20)

ROC = make_ROC_curve(outputs, 25) 
PR = make_precision_recall_curve(outputs, fontsize=25)


'''
One feature less -  No Probability
'''

all_features= ['IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein1', 'DomainFreqinProteome1']
all_features_renamed= ['IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']
exclude_feature= []

exclude_feature= ['Probability']

X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []
for boot in range(10): #10 random splits
    output_1= split_fit_rf(X_2, y_1, tree_n=100)
    outputs.append(output_1)

make_confusion_matrix(outputs)
make_ROC_curve(outputs, 25)
make_precision_recall_curve(outputs, fontsize=25)
make_feature_importance_plot(outputs, all_features_renamed, fontsize=20)

ROC = make_ROC_curve(outputs, 25)
PR = make_precision_recall_curve(outputs, fontsize=25)


'''
4 Features from motif and domains
'''

all_features= ['IUPredShort','metazoa_RLC', 'DomainFreqinProteome', 'DomainFreqbyProtein']
all_features_renamed= ['IUPredShort','metazoa_RLC', 'DomainFreqinProteome', 'DomainFreqbyProtein']
exclude_feature= [ 'Probability', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore']

# 

X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []
for boot in range(10):
    output_1= split_fit_rf(X_2, y_1, tree_n=100)
    outputs.append(output_1)

make_confusion_matrix(outputs)
make_ROC_curve(outputs, 25)
make_precision_recall_curve(outputs, fontsize=25)
make_feature_importance_plot(outputs, all_features_renamed, fontsize=20)

'''
Two features - IUPred and RLC
'''

all_features= ['IUPredShort','metazoa_RLC']
all_features_renamed= ['IUPredShort','metazoa_RLC']
exclude_feature= [ 'Probability', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqinProteome', 'DomainFreqbyProtein']


X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []

for boot in range(10): #10 random splits
    output_1= split_fit_rf(X_2, y_1, tree_n=100)
    outputs.append(output_1)

make_confusion_matrix(outputs)
make_ROC_curve(outputs, 25)
make_precision_recall_curve(outputs, fontsize=25)
make_feature_importance_plot(outputs, all_features_renamed, fontsize=20)

ROC = make_ROC_curve(outputs, 25)
PR = make_precision_recall_curve(outputs, fontsize=25)


'''
SINGLE features
'''

### IUPred
all_features= ['IUPredShort']
all_features_renamed= ['IUPredShort']
exclude_feature= [ 'Probability','metazoa_RLC', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqinProteome', 'DomainFreqbyProtein']

X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []
for boot in range(10): #10 random splits
    output_1= split_fit_rf(X_2, y_1, tree_n=100)
    outputs.append(output_1)

make_confusion_matrix(outputs)
make_ROC_curve(outputs, 25)
make_precision_recall_curve(outputs, fontsize=25)


AUC_d = []
AP_d = []

# Threshold test
for t in range(1, 101, 5):
    th =t/100
    X_2_thr = 1*(X_2>th)
    exclude_feature= []
    outputs = []
    for boot in range(10):
        output_1= split_fit_rf(X_2_thr, y_1, tree_n=100)
        outputs.append(output_1)
    
    AUC = make_ROC_curve(outputs, 25)
    AUC_d.append(AUC)
    AP = make_precision_recall_curve(outputs, fontsize=25)
    AP_d.append(AP)


### RLC
all_features= ['metazoa_RLC']
all_features_renamed= ['metazoa_RLC']
exclude_feature= [ 'Probability','IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqinProteome', 'DomainFreqbyProtein']


X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []
for boot in range(10):
    output_1= split_fit_rf(X_2, y_1, tree_n=100)
    outputs.append(output_1)

make_confusion_matrix(outputs)
make_ROC_curve(outputs, 25)
make_precision_recall_curve(outputs, fontsize=25)

ROC = make_ROC_curve(outputs, 25)
PR = make_precision_recall_curve(outputs, fontsize=25)


AUC_d = []
AP_d = []

# Threshold test 
for t in range(1, 101, 5):
    th =t/100
    X_2_thr = 1*(X_2>th)
    exclude_feature= []
    outputs = []
    for boot in range(10):
        output_1= split_fit_rf(X_2_thr, y_1, tree_n=100)
        outputs.append(output_1)
    
    AUC = make_ROC_curve(outputs, 25)
    AUC_d.append(AUC)
    AP = make_precision_recall_curve(outputs, fontsize=25)
    AP_d.append(AP)


'''
Logistic regression
'''
from sklearn.linear_model import LogisticRegression

all_features= ['IUPredShort','metazoa_RLC']
all_features_renamed= ['IUPredShort','metazoa_RLC']
exclude_feature= [ 'Probability', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqinProteome', 'DomainFreqbyProtein']

# 

X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []


rf, X_train, X_test, y_train, y_test= split_fit_rf(X_2, y_1)
# Logistic Regression model initialization
model = LogisticRegression(
    penalty='l1', #lasso
    C=2.0,
    solver='liblinear',
    max_iter=1000
)


# model training
model.fit(X_train, y_train)

# prédictions
y_pred = model.predict(X_test)

# model evaluation
accuracy = accuracy_score(y_test, y_pred)

print("Accuracy:", accuracy)

model.coef_[0]
X_2.columns

fig, ax= plt.subplots(figsize= (8,8))
roc = RocCurveDisplay.from_estimator(model, X_test, y_test, ax= ax)

fig, ax= plt.subplots(figsize= (8,8))
PrecisionRecallDisplay.from_estimator(model, X_test, y_test, ax= ax)

# Just IUPred
all_features= ['IUPredShort']
all_features_renamed= ['IUPredShort']
exclude_feature= [ 'Probability','metazoa_RLC', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqinProteome', 'DomainFreqbyProtein']

# 

X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []

rf, X_train, X_test, y_train, y_test= split_fit_rf(X_2, y_1)
# Logistic Regression model initialization
model = LogisticRegression(
    penalty='l1', #lasso
    C=2.0,
    solver='liblinear',
    max_iter=1000
)


# model training
model.fit(X_train, y_train)

model.coef_[0]
X_2.columns

fig, ax= plt.subplots(figsize= (8,8))
roc = RocCurveDisplay.from_estimator(model, X_test, y_test, ax= ax)

fig, ax= plt.subplots(figsize= (8,8))
PrecisionRecallDisplay.from_estimator(model, X_test, y_test, ax= ax)


# Just RLC
all_features= ['metazoa_RLC']
all_features_renamed= ['metazoa_RLC']
exclude_feature= [ 'Probability','IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqinProteome', 'DomainFreqbyProtein']

# 

X_2= X_1.drop(labels= exclude_feature, axis= 1)
exclude_feature= []
outputs = []

rf, X_train, X_test, y_train, y_test= split_fit_rf(X_2, y_1)
# Logistic Regression model initialization
model = LogisticRegression(
    penalty='l1', #lasso
    C=2.0,
    solver='liblinear',
    max_iter=1000
)


# model training
model.fit(X_train, y_train)

model.coef_[0]
X_2.columns

fig, ax= plt.subplots(figsize= (8,8))
roc = RocCurveDisplay.from_estimator(model, X_test, y_test, ax= ax)

fig, ax= plt.subplots(figsize= (8,8))
PrecisionRecallDisplay.from_estimator(model, X_test, y_test, ax= ax)

'''
Shared motif classes
'''

p_classes = PRS['Accession'].unique()
r_classes = RRS_1['Accession'].unique()

# Common motif classes
PRS_b = PRS[PRS.Accession.isin(r_classes)]
RRS_b = RRS_1[RRS_1.Accession.isin(p_classes)]

all_features= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein1', 'DomainFreqinProteome1']
all_features_renamed= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']
exclude_feature= []

df_1, X_1, y_1= preprocessing_dataset(PRS_b, RRS_b, all_features, all_features_renamed)

exclude_feature= []
outputs = []
for boot in range(10):  #10 random splits
    output_1= split_fit_rf(X_1, y_1)
    outputs.append(output_1)

make_confusion_matrix(outputs)
ROC_val = make_ROC_curve(outputs, 25)
PR_val = make_precision_recall_curve(outputs, fontsize=25)
make_feature_importance_plot(outputs, all_features_renamed, fontsize=20)

tprs= []
aucs= []
mean_fpr= np.linspace(0, 1, 100)
precisions= []
aps= []
mean_recall= np.linspace(0.0, 1.0, 100)

for n in range(10):  #Testins against to all 10 RF fits
    rf, X_train, X_test, y_train, y_test = outputs[n]
    fig, ax = plt.subplots(figsize= (8,8))
    prc = PrecisionRecallDisplay.from_estimator(rf, X_test, y_test, ax= ax)
    
    # Exclusive motif classes
    PRS_n = PRS[~ PRS.Accession.isin(r_classes)]
    RRS_n = RRS_1[~ RRS_1.Accession.isin(p_classes)]
    
    df_1, X_1, y_1= preprocessing_dataset(PRS_n, RRS_n, all_features, all_features_renamed)
    
    fig, ax = plt.subplots(figsize= (8,8))
    roc = RocCurveDisplay.from_estimator(rf, X_1, y_1, ax= ax)
    interp_tpr = np.interp(mean_fpr, roc.fpr, roc.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(roc.roc_auc)
    
    
    fig, ax = plt.subplots(figsize= (8,8))
    prc = PrecisionRecallDisplay.from_estimator(rf, X_1, y_1, ax= ax)
    interp_prec= np.interp(mean_recall, np.flipud(prc.recall), np.flipud(prc.precision))
    interp_prec[0] = 1.0
    interp_prec[-1] = len(y_test[y_test == 1])/ len(y_test)
    precisions.append(interp_prec)
    aps.append(prc.average_precision)

# Mean values for all extra-testing
mean_tpr = np.mean(tprs, axis= 0)
mean_tpr[-1] = 1.0
std_tpr = np.std(tprs, axis= 0)
mean_auc = np.mean(aucs)

mean_precision = np.mean(precisions, axis= 0)
std_precision = np.std(precisions, axis= 0)
mean_ap = np.mean(aps)
