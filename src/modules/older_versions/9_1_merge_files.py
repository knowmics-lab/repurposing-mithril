#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:09:44 2025

@ author: L-F-S
to run after connectivity score is done, to merge all files
"""
#%%
from conf import CS_OUT
import pandas as pd
import os
import subprocess
from preprocessing_utils import get_drugs_list



#%% Get all drugs
drugs_list=get_drugs_list()

# set mithril flag (0:DEG data, 1:mithril data)
mith=1

# set n of cores
cores=4

#%%
   

# get chunk size
chunk_size=int(len(drugs_list)/cores)
last_chunk_size=len(drugs_list)%cores

def get_chunk_indexes(i,chunk_size):
    i1=chunk_size*i
    i2=chunk_size*i+chunk_size
    return i1, i2
#%% merge files
dfs=[]
for i in range(cores):
    i1,i2=get_chunk_indexes(i, chunk_size)
    connectivity_dataset_filename=CS_OUT+str(i1)+'_'+str(i2)+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(i1)+'_'+str(i2)+'_mith_connectivity_score.tsv'

    dfs.append(pd.read_csv(connectivity_dataset_filename, sep='\t'))

connectivity_dataset_filename=CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_mith_connectivity_score.tsv'
dfs.append(pd.read_csv(connectivity_dataset_filename, sep='\t'))




connectivity_dataset_filename=CS_OUT+'DEG_connectivity_score.tsv' if mith==0 else CS_OUT+'mith_connectivity_score.tsv'

pd.concat(dfs).to_csv(connectivity_dataset_filename, sep='\t', index=False)
#%% delete files
for i in range(cores):
    i1,i2=get_chunk_indexes(i, chunk_size)
    connectivity_dataset_filename=CS_OUT+str(i1)+'_'+str(i2)+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(i1)+'_'+str(i2)+'_mith_connectivity_score.tsv'
    os.remove(connectivity_dataset_filename)

connectivity_dataset_filename=CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_DEG_connectivity_score.tsv' if mith==0 else CS_OUT+str(len(drugs_list)-chunk_size)+'_'+str(len(drugs_list))+'_mith_connectivity_score.tsv'
os.remove(connectivity_dataset_filename)