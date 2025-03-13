#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:41:12 2024

@ author: L-F-S
"""

#%%
import pandas as pd
from conf import HOME_DIR, BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG, TSR_OUT_DRUG

import subprocess
ex_mith=pd.read_csv(MITH_OUT_DRUG+'oligomycin-a_24h.perturbations.txt', sep='\t')
genes_list=ex_mith['Gene Id'].unique()
cores=30
print('n of genes:',len(genes_list), 'n of cores:', cores)

#%%
   
chunk_size=int(len(genes_list)/cores)
last_chunk_size=len(genes_list)%cores

def get_chunk_indexes(i,chunk_size):
    i1=chunk_size*i
    i2=chunk_size*i+chunk_size
    return i1, i2
#%%
for i in range(cores):
    i1,i2=get_chunk_indexes(i, chunk_size)
    subprocess.Popen(['python','map_mith3_drug_wise_output_to_gene_wise_metanalysis.py',str(i1), str(i2)])

subprocess.Popen(['python','map_mith3_drug_wise_output_to_gene_wise_metanalysis.py',str(len(genes_list)-chunk_size),str(len(genes_list))])
