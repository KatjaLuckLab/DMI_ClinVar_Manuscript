#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 16:22:25 2025

"""


'''
Lbraries
'''

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



'''
Functions adapted from Chop Yan Lee DMI Predictor

'''
def preprocessing_dataset(PRS_input, RRS_input, all_features, all_features_renamed): 
    """
    Takes the PRS and RRS, concatenate them and preprocessing the NaNs and dummy value.
    
    Args:
        PRS_path (str): Absolute path to the PRS dataset
        RRS_path (str): Absolute path to the RRS dataset

    Returns:
        df (pd.DataFrame): The concatenated PRS and RRS without the rows with NaNs and dummy value
        X (pd.DataFrame): The features
        y (pd.DataFrame): The outcomes
    """
    
    PRS_input['label']= 1
    RRS_input['label']= 0
    for df in [PRS_input, RRS_input]:
        df.replace(88888, df.DomainEnrichment_zscore.median(), inplace= True)
        for ind, row in df.iterrows():
            if pd.notna(row['DomainFreqbyProtein2']):
                df.loc[ind, 'DomainFreqbyProtein1'] = np.mean([row['DomainFreqbyProtein1'], row['DomainFreqbyProtein2']])
                df.loc[ind, 'DomainFreqinProteome1'] = np.mean([row['DomainFreqinProteome1'], row['DomainFreqinProteome2']])
    df= pd.concat([PRS_input, RRS_input], axis= 0, ignore_index= True)
    df.dropna(subset= all_features, inplace= True)
    df.rename(columns= {'DomainFreqbyProtein1': 'DomainFreqbyProtein', 'DomainFreqinProteome1': 'DomainFreqinProteome'}, inplace= True)
    X= df[all_features_renamed].copy()
    y= df['label']
    return df, X, y 





'''
Random forest
'''
def split_fit_rf(X, y, exclude_feature= None, tree_n =100): # features to be dropped saved as list
    """
    Split the provided features and outcomes into train and test set and fit a random forest classifier to the train set.

    Args:
        X (pd.DataFrame): The features
        y (pd.DataFrame): The outcomes
        exclude_feature (list of str): Features to be excluded, default None

    Returns:
        rf (sklearn.ensemble.RandomForestClassifier): The fitted random forest model
        X_train (pd.DataFrame): The feature values in the train set
        X_test (pd.DataFrame): The feature values in the test set
        y_train (pd.DataFrame): The outcome in the train set
        y_test (pd.DataFrame): The outcome in the test set
    """
    rf= RandomForestClassifier(n_estimators= tree_n, oob_score= True, verbose= True, n_jobs= 1) # Tested the model again on 11.03.2024 and found that it only works if n_jobs is set to 1. -1 makes use of all CPU cores but somehow it does not work on my computer.

    if exclude_feature != None:
        X= X.drop(labels= exclude_feature, axis= 1)
    X_train, X_test, y_train, y_test= train_test_split(X, y,  stratify= y, test_size=0.2) #random_state= 0,
    rf.fit(X_train, y_train)

    return rf, X_train, X_test, y_train, y_test



'''
Evaluation
'''
def make_confusion_matrix(split_fit_outputs):
    """
    Plot the confusion matrix of the fitted model applied on the test set using the returned variables from the function split_fit_rf and save the figure as pdf file.

    Args:
        split_fit_outputs (list): List of three outputs (that correspond to triplicates of each RRS version) from the function split_fit_rf
    """
    fig, axes= plt.subplots(1, 10, figsize= (10, 10))

    for i, ele in enumerate(zip(split_fit_outputs, axes)):
        output, ax= ele
        rf, X_train, X_test, y_train, y_test= output
        ConfusionMatrixDisplay.from_estimator(rf, X_test, y_test, ax= ax, colorbar= False)
        ax.set_title(f'RRSv1_{i+1}', fontsize = 10)



def make_ROC_curve(split_fit_outputs, fontsize= 12):
    """
    Plot the individual ROC curves in one figure, as well as the averaged ROC curves across triplicates +/- 1 std. in another figure, of the fitted model applied on the test set using the returned variables from the function split_fit_rf and save the figures as pdf file. The average TPR, FPR and their standard deviations are written out in a tsv file.

    Args:
        split_fit_outputs (list): List of three outputs (that correspond to triplicates of each RRS version) from the function split_fit_rf
    """
    tprs= []
    aucs= []
    mean_fpr= np.linspace(0, 1, 100)

    fig, ax= plt.subplots(figsize= (8,8))

    for i, output in enumerate(split_fit_outputs):
        rf, X_train, X_test, y_train, y_test= output
        roc = RocCurveDisplay.from_estimator(rf, X_test, y_test, ax= ax)
        interp_tpr = np.interp(mean_fpr, roc.fpr, roc.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(roc.roc_auc)

    # plot avg ROC curve across triplicates
    mean_tpr= np.mean(tprs, axis= 0)
    mean_tpr[-1]= 1.0
    std_tpr= np.std(tprs, axis= 0)
    mean_auc= np.mean(aucs)
    

    fig, ax= plt.subplots(figsize= (8,8))

    ax.plot(mean_fpr, mean_tpr, label= f'mean AUC = {mean_auc:.2f}')
    ax.fill_between(mean_fpr, mean_tpr - std_tpr, mean_tpr + std_tpr, alpha= 0.3, label= '1 std. dev.')
    ax.set(xlim=[-0.05, 1.05], ylim= [-0.05, 1.05])
    ax.legend(loc= 'lower right', fontsize = fontsize*0.8)
    ax.tick_params(axis='both', which='major', labelsize=fontsize*0.6)
    ax.set_xlabel('False Positive Rate (Positive label: 1)', fontsize = fontsize*0.8)
    ax.set_ylabel('True Positive Rate (Positive label: 1)', fontsize = fontsize*0.8)
    ax.set_title(f'ROC curve averaged across {len(split_fit_outputs)} samples', fontdict= {'fontsize': fontsize})
    
    return(mean_fpr, mean_tpr)


def make_precision_recall_curve(split_fit_outputs, fontsize=12):
    """
    Plot the individual PR curves in one figure, as well as the averaged PR curves across triplicates +/- 1 std. in another figure, of the fitted model applied on the test set using the returned variables from the function split_fit_rf and save the figures as pdf file. The average recall, precision and their standard deviations are written out in a tsv file.

    Args:
        split_fit_outputs (list): List of three outputs (that correspond to triplicates of each RRS version) from the function split_fit_rf
    """
    precisions= []
    aps= []
    mean_recall= np.linspace(0.0, 1.0, 100)

    fig, ax= plt.subplots(figsize= (8,8))

    for i, output in enumerate(split_fit_outputs):
        rf, X_train, X_test, y_train, y_test= output
        prc = PrecisionRecallDisplay.from_estimator(rf, X_test, y_test, ax= ax)
        interp_prec= np.interp(mean_recall, np.flipud(prc.recall), np.flipud(prc.precision))
        interp_prec[0]= 1.0
        interp_prec[-1]= len(y_test[y_test == 1])/ len(y_test)
        precisions.append(interp_prec)
        aps.append(prc.average_precision)

    # plot avg PR curve across triplicates
    mean_precision= np.mean(precisions, axis= 0)
    # mean_precision[-1]= 0.5
    std_precision= np.std(precisions, axis= 0)
    mean_ap= np.mean(aps)

    fig, ax= plt.subplots(figsize= (8,8))

    ax.plot(mean_recall, mean_precision, label= f'mean AP = {mean_ap:.2f}')
    ax.fill_between(mean_recall, mean_precision - std_precision, mean_precision + std_precision, alpha= 0.3, label= '1 std. dev.')
    ax.set(xlim=[-0.05, 1.05], ylim= [0.45, 1.05])
    ax.legend(loc= 'lower left', fontsize = fontsize*0.8)
    ax.tick_params(axis='both', which='major', labelsize=fontsize*0.6)
    ax.set_xlabel('Recall (Positive label: 1)', fontsize= fontsize*0.8)
    ax.set_ylabel('Precision (Positive label: 1)', fontsize= fontsize*0.8)
    ax.set_title(f'PR curve averaged across {len(split_fit_outputs)} samples', fontdict= {'fontsize': fontsize})
    
    return(mean_recall, mean_precision)



def make_cvacc_oob_acc_plot(split_fit_outputs, fontsize= 12):
    """
    Plot a barplot showing the accuracy of models fitted to the individual replicate of an RRS version. The accuracies of models are evaluated using three strategies, 5 fold cross validation accuracy, out-of-bag samples, and accuracy evaluated on the test set. The figure is saved as pdf file.

    Args:
        split_fit_outputs (list): List of three outputs (that correspond to triplicates of each RRS version) from the function split_fit_rf
    """
    cvacc= []
    cvacc_std= []
    oob= []
    acc= []

    for i, output in enumerate(split_fit_outputs):
        rf, X_train, X_test, y_train, y_test= output
        oob.append(rf.oob_score_)
        acc.append(accuracy_score(y_test, rf.predict(X_test)))
        cvacc.append(np.mean(cross_val_score(RandomForestClassifier(n_estimators= 1000), X_train, y_train, cv= 5, n_jobs= -1)))
        cvacc_std.append(np.std(cross_val_score(RandomForestClassifier(n_estimators= 1000), X_train, y_train, cv= 5, n_jobs= -1)))

    N= 3
    ind= np.arange(N)
    width= 0.25

    plt.figure(figsize= (8,6))

    plt.bar(ind, cvacc, width, yerr= cvacc_std, capsize= 7, color= 'c', label= '5-fold CV')
    plt.bar(ind + width, oob, width, color= 'b', label= 'oob')
    plt.bar(ind + 2 * width, acc, width, color= 'k', label= 'holdout')

    plt.xticks(ind + width, [f'r_1', f'r_2', f'r_3'], fontsize= fontsize)
    plt.ylim(0, 1)
    plt.title(f'Comparison of accuracy scores computed with different strategies ', fontsize= fontsize)
    plt.ylabel('Accuracy score', fontsize= fontsize)
    plt.legend(bbox_to_anchor= (1.05, 1), loc= 'upper left', borderaxespad= 0.)
    # plt.grid(alpha= 0.2)





def make_feature_importance_plot(split_fit_outputs, all_features_renamed, exclude_feature= None, fontsize=12):
    """
    Plot a horizontal barplot showing the feature importance of models fitted to the individual replicate of an RRS version and another horizontal barplot showing the average feature importance of fitted models to all replicates of an RRS version. The feature importance is evaluated using Gini index. The figures are saved as pdf files.

    Args:
        split_fit_outputs (list): List of three outputs (that correspond to triplicates of each RRS version) from the function split_fit_rf
        exclude_feature (list of str): Features to be excluded, default None
    """
    for i, output in enumerate(split_fit_outputs):
        rf, _, _, _, _= output
        if i== 0:
            if exclude_feature != None:
                features= list(filter(lambda feat: feat not in exclude_feature, all_features_renamed))
                xlim= None
            else:
                features= all_features_renamed
                xlim= [0, 0.5]
            feat_imp_df= pd.DataFrame(data= {'Features': features, f'R_{i+1}': rf.feature_importances_})
        else:
            feat_imp_df= pd.concat([feat_imp_df, pd.Series(rf.feature_importances_, name= f'R_{i+1}')], axis= 1)
    
    sets = len(split_fit_outputs)
    max_feature = feat_imp_df.iloc[:, 1:sets+1].max().max()
    feat_imp_df['mean']= feat_imp_df.iloc[:, 1:sets+1].mean(axis= 1, numeric_only=True)
    feat_imp_df['std']= feat_imp_df.iloc[:, 1:sets+1].std(axis= 1, numeric_only=True)
    feat_imp_df= feat_imp_df.sort_values(by= 'R_1', ascending= True)

    N= len(features)
    ind= np.arange(N)
    height= 0.3

    plt.figure(figsize= (8,10))
    for rrs in feat_imp_df:
        plt.barh(ind + 2*height, feat_imp_df[rrs], height, color= 'c', label= 'rrs_1')
        plt.barh(ind + height, feat_imp_df[rrs], height, color= 'b', label= 'rrs_2')
        plt.barh(ind, feat_imp_df[rrs], height, color= 'r', label= 'rrs_3')

    plt.yticks(ind + height, feat_imp_df.Features)
    plt.legend(loc= 'best')
    plt.title(f'Feature importance of rrs')
    plt.xlabel('Gini index')
    # plt.grid(alpha= 0.2)


    # plot avg and std of feature importance across RRS triplicates.
    feat_imp_df= feat_imp_df.sort_values(by= 'mean', ascending= True)

    plt.figure(figsize= (6, 8))

    plt.barh(feat_imp_df['Features'], feat_imp_df['mean'], alpha= 0.9, color= 'steelblue', xerr= feat_imp_df['std'], capsize= 5)

    plt.xlim([0,max_feature*1.1])
    plt.title(f'Avg of feature importance across {len(split_fit_outputs)} samples ', fontsize= fontsize)
    plt.xlabel('Gini index', fontsize= fontsize*0.8)    
    plt.tick_params(axis='both', which='major', labelsize=fontsize*0.6)
    # plt.grid(alpha= 0.2)

