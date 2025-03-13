#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 09:50:20 2024

@ author: L-F-S

Maps row names of genes of Mithril 3 input matrix from gene names to gene ids

Creates dictionary of gene symbols to unique gene ids and
saves DEG disease signatures with unique gene ids

TO do the same for DEG drug signatures, run the apporpriate cell 
in preproc_utils.py. Not the most straight forward way, but it works
"""

import os
import sys
import numpy as np
import pandas as pd
import pickle
import collections
from conf import BASE_DIR, MITH_IN_DRUG, MITH_OUT_DRUG, DISEASE,\
MITH_IN_DISEASE, TSR_OUT_DISEASE, symbol_2geneid, alias_2geneid

#%% DRUG DATA
LINCS_metanalysis_filename =MITH_IN_DRUG+'LINCS_metanalysis_genename.mi'
LINCS_metanalysis_matrix=pd.read_csv(LINCS_metanalysis_filename,sep='\t',header=0)
#%%

def map_name_to_id(genename):
    if genename in symbol_2geneid.keys():
        return symbol_2geneid[genename]
    return genename
LINCS_metanalysis_matrix['Unnamed: 0']=LINCS_metanalysis_matrix['Unnamed: 0'].apply(lambda x : map_name_to_id(x))
    #%%
LINCS_metanalysis_filename_id='LINCS_metanalysis.mi'
LINCS_metanalysis_matrix.to_csv(MITH_IN_DRUG+LINCS_metanalysis_filename_id,sep='\t', index=False)

#%%

######################
# DISEASE DATA
#####################
# preprocessing to remove duplicates

disease_gene_symbol_mith_in=pd.read_csv(MITH_IN_DISEASE+DISEASE+'_signature.mi', header=None, sep='\t', decimal=',')
disease_genes=disease_gene_symbol_mith_in[0]
#%%
genes_not_in_alias=[]
disease_genes_id_to_symbol_univocous={}
duplicated_ids_with_different_symbols=collections.defaultdict(list)
for symbol in disease_genes:
    try:
        gene_id=str(alias_2geneid[symbol])
#        if not gene_id in disease_genes_id_to_symbol_univocous.keys():
        disease_genes_id_to_symbol_univocous[gene_id] = symbol #actually not univocous yet, there are duplicates
        duplicated_ids_with_different_symbols[gene_id].append(symbol) # 480 GENI DUPLICATI nell alias
    except:

        genes_not_in_alias.append(symbol)
print(len(genes_not_in_alias), 'gene symbols without gene id ') #198 genes not in alias

#only keep duplicates:
topop=[]
for gene_id, symbols  in duplicated_ids_with_different_symbols.items():
    if len(symbols)==1:
        topop.append(gene_id)
for gene_id in topop:
    duplicated_ids_with_different_symbols.pop(gene_id)
print(len(duplicated_ids_with_different_symbols.keys()), 'duplicated ids')

#%%

def print_duplicated(genes,df):
    FCs=[]
    for gene_symbol in genes:
        fc=disease_gene_symbol_mith_in[disease_gene_symbol_mith_in[0]==gene_symbol][1].iloc[0]
#        np.round(fc,3)
#        print(fc, type(fc))
        FCs.append(np.round(fc,3))
    if not len(np.unique(np.array(FCs)))==1:
#        print('duplicated:',genes)

        return 1
    return 0

count_duplicated_genes_with_different_FC=0
for gene_id, duped_symbols in duplicated_ids_with_different_symbols.items():
    count_duplicated_genes_with_different_FC+=print_duplicated(duped_symbols, disease_gene_symbol_mith_in)
print('tot duplicated genes with idfferent FC', count_duplicated_genes_with_different_FC)
#il problema e che hanno dati diversi, ma noi faremo finta di niente
# e li leviamo proprio di torno questi geni.
# topop è la lista di geni univoci sono 15k
genes_id_to_symbol_univocous={gene_id:disease_genes_id_to_symbol_univocous[gene_id] for gene_id in topop}

#%%save dictionary for later translation:
filename=BASE_DIR+'other_data/'+DISEASE+'_gene_id_to_symbol.pkl'
with open(filename, 'wb') as f:
    pickle.dump(genes_id_to_symbol_univocous, f)
    #%% DONE AFTER THIS PART, WE CAN TRANSLATE THE MITHRIL INPUT, by keeping only relevant genes
gene_id_mith_input=[]
for gene_id in topop: #lentino, ma tant è lo faccio una volt asoloa
    fc=disease_gene_symbol_mith_in[disease_gene_symbol_mith_in[0]==disease_genes_id_to_symbol_univocous[gene_id]][1].iloc[0]
    gene_id_mith_input.append(gene_id+'\t'+str(fc))

#%% save mithril input
    
f=open(MITH_IN_DISEASE+DISEASE+'_signature_gene_id.mi','w')
f.write(('\n').join(gene_id_mith_input))
f.close()

#%%% with the chosen dictionary, convert the original DEG data
# not the most straightforward way, it would have been more straightforward
# to fist convert the DEG input once, but hey, it works.

filename=BASE_DIR+'other_data/'+DISEASE+'_gene_id_to_symbol.pkl'
with open(filename, 'rb') as f:
    genes_id_to_symbol_univocous=pickle.load(f)

DEG_disease_signature= pd.read_csv(TSR_OUT_DISEASE+DISEASE+'/'+DISEASE+'_signature.csv', sep=';',decimal=',')

# remove gene symbols not in mapping
DEG_disease_signature = DEG_disease_signature[~DEG_disease_signature.gene.isin(genes_not_in_alias)]

# Remove duplicates
DEG_disease_signature['gene_id']=DEG_disease_signature['gene'].apply(lambda gene : str(alias_2geneid[gene]))
DEG_disease_signature = DEG_disease_signature[DEG_disease_signature.gene_id.isin(topop)]

# Save file
DEG_disease_signature.to_csv(TSR_OUT_DISEASE+DISEASE+'/'+DISEASE+'_signature_gene_id.csv', sep=';',decimal=',', index=None)
