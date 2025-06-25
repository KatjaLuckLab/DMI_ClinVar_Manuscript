#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:25:09 2024

@author: jesusav
"""

import os
import numpy as np
import re
import argparse
import pandas as pd
import warnings


'''
Make a list of all proteins in the proteome
IDs are the dictionary keys
'''

os.chdir('/Volumes/imb-luckgr/imb-luckgr2/projects/AlphaFold/HuRI_DMI_AF_modelling/Domain_mappings/')
#os.chdir('.')

HumanProteome = 'DMI_seqs.fasta' #File with information on Proteins 
    
proteins_file = open(HumanProteome, 'r') #open up file
proteins_seqs = {}

for line in proteins_file: #for every line in the protein file     
    if line.startswith(">"):
        row = line.strip()
        ID = row[1:] #get the ID
        
        proteins_seqs[ID] = [""]
    else:
        row = line.strip()
        proteins_seqs[ID][0] = proteins_seqs[ID][0] + row #The name skiping the column of SigPep prob
        
del row, line, ID


'''
Get DMI boundary information and select unique protein pair ids
'''

os.chdir('/Volumes/imb-luckgr/imb-luckgr2/projects/AlphaFold/HuRI_DMI_AF_modelling/Domain_mappings/')
#os.chdir('.')

# DMI_boundaries
DMI_boundaries = pd.read_table('HuRI_human_DMI_Clinvar_PU_highScr_noMOD_repeated.tsv')

DMI_boundaries = DMI_boundaries.dropna(subset=['pLDDT_start'])
DMI_proteins = DMI_boundaries["DomainProtein"].unique().tolist() + DMI_boundaries["SLiMProtein"].unique().tolist()
  

'''
Produce table with fasta files limits
'''
out_file = open('DMI_fragment_boundary_entensions_repeated.txt',"w")
out_file.write('DomainProtein\tDomainID1\tdomain_start\tdomain_end\tmotifID\tAccession\tELM_start\tELM_end\tmotif_start\tmotif_end\n' ) 

for index, row in DMI_boundaries.iterrows():
    domainProtein = row.DomainProtein # Domain ID
    domain_seq = proteins_seqs[domainProtein][0][int(row.pLDDT_start-1):int(row.pLDDT_end-1)] # Domain sequence
    domain_start = int(row.pLDDT_start)  
    domain_end = int(row.pLDDT_end) 
    
    motifID = row.SLiMProtein # Motif ID
    motif_start = row.ELM_start - 1 # Motifs do not have python index defined boundaries, so to get their sequence we need -1
    motif_end = row.ELM_end - 1 
    if motif_start > 4:
        m_s = int(motif_start)-5
    else:
        m_s = int(motif_start)
        
    motif_seq = proteins_seqs[motifID][0][m_s:(int(motif_end)+5)] # Motif sequence plus flanking regions
    
    out_file.write(f'{domainProtein}\t{row.DomainID1}\t{domain_start}\t{domain_end}\t{motifID}\t{row.Accession}\t{row.ELM_start}\t{row.ELM_end}\t{m_s}\t{int(motif_end)+5}\n' ) 

out_file.close()


'''
Produce fasta files with protein fragment pairs
'''

# os.chdir('/Volumes/imb-luckgr/imb-luckgr2/projects/AlphaFold/HuRI_DMI_AF_modelling/fastaFiles')
os.chdir('/Volumes/imb-luckgr/imb-luckgr2/projects/AlphaFold/HuRI_DMI_AF_modelling/missing_fastas')

for index, row in DMI_boundaries.iterrows():
    domainProtein = row.DomainProtein # Domain ID
    domain_seq = proteins_seqs[domainProtein][0][int(row.pLDDT_start-1):int(row.pLDDT_end-1)] # Domain sequence
    domain_start = int(row.pLDDT_start)  
    domain_end = int(row.pLDDT_end) 
    
    motifID = row.SLiMProtein # Motif ID
    motif_start = row.ELM_start - 1 # Motifs do not have python index defined boundaries, so to get their sequence we need -1
    motif_end = row.ELM_end - 1 
    if motif_start > 4:
        m_s = int(motif_start)-5
    else:
        m_s = int(motif_start)
        
    motif_seq = proteins_seqs[motifID][0][m_s:(int(motif_end)+6)] # Motif sequence plus flanking regions
    
    out_file = open(f'{domainProtein}_{row.DomainID1}_O_{domain_start}_{domain_end}.{motifID}_{row.Accession}_D_{row.ELM_start}_{row.ELM_end}.fasta',"w")
    out_file.write(f'>{domainProtein}_{row.DomainID1}_O_{domain_start}_{domain_end}\n' ) 
    out_file.write(f'{domain_seq}\n' ) 
    out_file.write(f'>{motifID}_{row.Accession}_D_{row.ELM_start}_{row.ELM_end}\n' ) 
    out_file.write(f'{motif_seq}' ) 
    out_file.close()






