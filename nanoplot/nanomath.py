# wdecoster
import numpy as np
import pandas as pd


def getN50(a):
	'''
	Calculator function: Get read N50.
	Based on https://github.com/PapenfussLab/Mungo/blob/master/bin/fasta_stats.py
	'''
	return a[np.where(np.cumsum(a)>=0.5*np.sum(a))[0][0]]


def removeLengthOutliers(df, columnname):
	'''
    Calculation function: Remove records with length-outliers
    '''
	return df[df[columnname] < (np.median(df[columnname]) + 3 * np.std(df[columnname]))]


def aveQual(quals):
	'''	Calculation function: Receive the integer quality scores of a read and return the average quality for that read'''
	return sum(quals) / len(quals)
