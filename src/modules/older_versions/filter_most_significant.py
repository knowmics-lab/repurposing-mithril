#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 18:45:21 2025

@ author: L-F-S
"""
import os
if not os.getcwd().endswith('modules'):
    os.chdir('modules')
import sys
import numpy as np
import pandas as pd
import pickle
from conf import HOME_DIR, BASE_DIR, DISEASE,MITH_IN_DRUG, MITH_OUT_DRUG,\
 TSR_OUT_DRUG, TSR_OUT_DISEASE, TSR_OUT_CSCORE, alias_2geneid,load_disease_signature
 
print(DISEASE)
 #%%
mith_filename='als_NYGC_mithril_signature.csv'
filename=TSR_OUT_DISEASE+DISEASE+'/'+mith_filename
disease_signature=pd.read_csv(filename, sep=';',decimal=',')
#%%

#load various significance files
150_most_significant_genes_filename <- paste0(als_NYGC_cfg$signature_most_significant_genes_base_dir, "als_NYGC_150_most_significant_genes.csv")
als_NYGC_cfg$signature_150_most_significant_landmark_genes_filename <- paste0(als_NYGC_cfg$signature_most_significant_genes_base_dir, "als_NYGC_150_most_significant_landmark_genes.csv")
als_NYGC_cfg$signature_bin_chen_most_significant_genes_filename <- paste0(als_NYGC_cfg$signature_most_significant_genes_base_dir, "als_NYGC_bin_chen_most_significant_genes.csv")


als_NYGC_cfg$connectivity_score_base_dir <- "output/connectivity_score/als_NYGC_lorenzo/"
als_NYGC_cfg$connectivity_score_150_MS_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_150_MS_genes.csv")
als_NYGC_cfg$connectivity_score_150_MS_landmark_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_150_MS_landmark_genes.csv")
als_NYGC_cfg$connectivity_score_bin_chen_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_bin_chen_genes.csv")

als_NYGC_cfg$connectivity_score_PGx_150_MS_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_PGx_150_MS_genes.csv")
als_NYGC_cfg$connectivity_score_PGx_bin_chen_genes_filename <- paste0(als_NYGC_cfg$connectivity_score_base_dir, "als_NYGC_connectivity_score_PGx_bin_chen_genes.csv")

