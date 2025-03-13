#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:10:47 2024

@ author: L-F-S
"""
#%%
import os
import sys
import numpy as np
import pandas as pd
from conf import MITH_IN_DRUG, MITH_OUT_DRUG

df_list=[]
for mi2 in os.listdir(MITH_IN_DRUG):
    drug_condition_name=mi2.rstrip('.mi')
    if not drug_condition_name.startswith('LINCS'):
        print(drug_condition_name)
        df=pd.read_csv(MITH_IN_DRUG+mi2,sep='\t',header=None,index_col=0, names=[drug_condition_name])
        print(df.shape)
        print(df.columns)
        df_list.append(df)

matrix=pd.concat(df_list, axis=1)
print(matrix.shape)
print(matrix.columns)
print(matrix.head(3))
matrix.to_csv(MITH_IN_DRUG+'LINCS_metanalysis_genename.mi', sep='\t')
