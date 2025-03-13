#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:00:15 2025

@ author: L-F-S
"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import pickle
from conf import HOME_DIR, BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG, TSR_OUT_DRUG, TSR_OUT_DRUG_META

filenames=[name for name in os.listdir(TSR_OUT_DRUG_META+'gene_id/') if name.endswith('.txt')]

#%%


for filename in filenames:
    print(filename)
    data=pd.read_csv(TSR_OUT_DRUG_META+'gene_id/'+filename, sep='\t', index_col=0)
    
    # get gene name from file
    genename=filename.split('_')[0]
    data['gene']=genename
    data.to_csv(TSR_OUT_DRUG_META+filename, sep='\t')