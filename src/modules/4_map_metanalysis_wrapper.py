#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:41:12 2024

@ author: L-F-S
"""

import subprocess



nfinale=1
for i in range(30):
    string='Rscript map_metanalysis_gene_wise_to_drug_wise.R ' + str(nfinale) + ' ' + str(nfinale+99)
    print(string)
    subprocess.Popen(['Rscript','map_metanalysis_gene_wise_to_drug_wise.R',str(nfinale), str(nfinale+99)])
    nfinale+=100
    

subprocess.Popen(['Rscript','map_metanalysis_gene_wise_to_drug_wise.R','3001','3222'])
