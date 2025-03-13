#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

Maps row names of genes of Mithril 3 input matrix from gene names to gene ids
"""
import os
import sys
import numpy as np
import pandas as pd
import pickle
from conf import HOME_DIR, BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG

LINCS_metanalysis_filename =MITH_IN_DRUG+'LINCS_metanalysis_genename.mi'
LINCS_metanalysis_matrix=pd.read_csv(LINCS_metanalysis_filename,sep='\t',header=0)
#%%
alias_2geneid_filename=BASE_DIR+'utils/alias_2geneid.pkl'
with open(alias_2geneid_filename, 'rb') as f:
    alias_2geneid = pickle.load(f)

id_2alias={y:x for (x,y) in alias_2geneid.items()}

def map_id_to_name(idd):
    if idd in id_2alias.keys():
        return id_2alias[idd]
    raise ValueError('No corresponding gene name for id', idd)
    return idd

LINCS_metanalysis_matrix['Unnamed: 0']=LINCS_metanalysis_matrix['Unnamed: 0'].apply(lambda x : map_id_to_name(x))
#%%
LINCS_metanalysis_filename_id='LINCS_metanalysis.mi'
LINCS_metanalysis_matrix.to_csv(MITH_IN_DRUG+LINCS_metanalysis_filename_id,sep='\t', index=False)