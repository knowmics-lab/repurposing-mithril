#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:41:12 2024

@ author: L-F-S
"""

#%%
import pandas as pd
from conf import HOME_DIR, BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG, TSR_OUT_DRUG
import os
import numpy as np
import subprocess
#%%
# Open a random file to get all mithril genes (and other things) list

ex_mith=pd.read_csv(MITH_OUT_DRUG+'oligomycin-a_24h.perturbations.txt', sep='\t')

print(ex_mith.columns)
#
##'# Pathway Id', 'Pathway Name', 'Gene Id', 'Gene Name', 'Perturbation', 'Accumulator', 'pValue'
print(ex_mith.shape)

# Remove special cahrs from gene name sto avoid erros in later file naming
def remove_special_characters(df):
    special_chars=';/:*?\"|'
    for char in special_chars:
        indexes=df[df['Gene Name'].str.contains('\\'+char)].index
        if len(indexes)>0:
            print(len(indexes),char)
            df.drop(index=indexes, inplace=True)
    return df

ex_mith=remove_special_characters(ex_mith)

#%% Printing ad file exploring
print(len(ex_mith['Gene Id'].unique()), 'genes')
# genes appear several times , let's check that they always have the same perturbation value


# printing  stuff to see wat's inside ex_mith
#for gene in ex_mith['Gene Id'].unique():
#    print(gene, len(ex_mith[ex_mith['Gene Id']==gene]), 'unique pert', \
#          len(ex_mith[ex_mith['Gene Id']==gene]['Perturbation'].unique()),\
#          'unique pval', len(ex_mith[ex_mith['Gene Id']==gene]['pValue'].unique()),\
#          'unique acc', len(ex_mith[ex_mith['Gene Id']==gene]['Accumulator'].unique()),\
#          'unique path', len(ex_mith[ex_mith['Gene Id']==gene]['# Pathway Id'].unique()))
#    break

# so this just prints the PATHWAYS in which a given gene appears!
          # important 8write in log!!|!
# every line is A GENE IN A PATHWAY, but its value of acc,pval and pert is the same in every pathway.
          # but other question
genes_list=np.sort(list(ex_mith['Gene Id'].unique()))
cores=30
print('n of genes:',len(genes_list), 'n of cores:', cores)

#%%
   
#genes already done
done_genes=[x.split('_')[1] for x in os.listdir(TSR_OUT_DRUG+'LINCS_lorenzo/metanalysis_mith3_gene_wise/')]


# get chunk size
chunk_size=int(len(genes_list)/cores)
last_chunk_size=len(genes_list)%cores

def get_chunk_indexes(i,chunk_size):
    i1=chunk_size*i
    i2=chunk_size*i+chunk_size
    return i1, i2

#%% 
for i in range(cores):
    i1,i2=get_chunk_indexes(i, chunk_size)
    subprocess.Popen(['python','map_mith3_drug_wise_output_to_gene_wise_metanalysis_iterative.py',str(i1), str(i2)])

subprocess.Popen(['python','map_mith3_drug_wise_output_to_gene_wise_metanalysis_iterative.py',str(len(genes_list)-chunk_size),str(len(genes_list))])
