#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:25:09 2024

@author: jesusav
"""

import os
import numpy as np
import re
import pandas as pd
import warnings

## This file takes a text file containing all AF2 pdb files on a directory and extract their pLDDT values
## it saves the values in two text files: the first with one residue per row and the second with a list of values per protein



'''
Make a list of all proteins in the proteome
IDs are the dictionary keys
Information store include gene name, seq length and sequence
'''

os.chdir('/Volumes/imb-luckgr/imb-luckgr2/projects/AlphaFold/AFDB_HumanReferenceProteome/')
#os.chdir('.')

PFAM_domains = pd.read_csv('HuRI_human_DMI_domains.csv')


PFAM_domains = pd.read_csv('HuRI_human_DMI_domains_repeated.csv')

# PFAM_domains = pd.read_table('high_pLDDTregions/HumanReferenceProteome_high_pLDDTregions_70.tsv')
# PFAM_domains = PFAM_domains.rename(columns={'Region_protein': 'DomainProtein', 'Region_n': 'DomainID1', 'region_start': 'domain_start', 'region_end': 'domain_end'})

PFAM_domains.insert(0, "index", range(0, len(PFAM_domains), 1), True)
PFAM_domains_list = PFAM_domains.set_index('index').transpose().to_dict('list') #put the datafram in dictionary (of lists) format



plddt_table = pd.read_table('proteins_pLDDT_values_list.txt')
plddt_table[['UniProt', 'Fragment']] = plddt_table['UniProt_ID'].str.split('-', expand=True)
plddt_table = plddt_table[(plddt_table['Fragment'] == "F1")]
plddt_table = plddt_table.drop(['UniProt_ID', 'Fragment'], axis="columns")

plddt_list_raw = plddt_table.set_index('UniProt').transpose().to_dict('list') #put the datafram in dictionary (of lists) format
plddt_list = {}
for item in plddt_list_raw:
    temp = plddt_list_raw[item][0]
    temp = re.sub("\[|\]", "", temp)
    temp = temp.split(",")
    temp = [ float(n) for n in temp ]
    plddt_list[item] = temp
del item, temp, plddt_list_raw




'''
Extension parameters
'''
border_window = 10 # Residue window in which the pLDDT values would be averaged, checked and compared
pLDDT_border_thrs = 80 # pLDDT value that the domain window has to overcome in order to be extended (to be kept)


'''
Extension functions
'''

### AlphaFold has individual fragment pdb files

def check_boundary(pLDDT_VAL, DOM_B, B_TYPE): # B - boundary
    DOM_B = DOM_B - 1    
    if B_TYPE=="Start":
        with warnings.catch_warnings():
          warnings.simplefilter("ignore", category=RuntimeWarning) 
          DOM_S_I_AVG = np.nanmean(pLDDT_VAL[DOM_B:(DOM_B+border_window)]) # I - Inner, S - Start using nanmean as there might be empty values in the list if values
         
        if DOM_S_I_AVG < pLDDT_border_thrs: # wheter the start inner window
            return True
        else:
            return False
        
    elif B_TYPE=="End":
        with warnings.catch_warnings():
          warnings.simplefilter("ignore", category=RuntimeWarning)
          DOM_E_I_AVG = np.nanmean(pLDDT_VAL[(DOM_B-border_window):DOM_B]) # E - End, I - Inner
                    
        if DOM_E_I_AVG < pLDDT_border_thrs: # wheter the start inner window
            return True
        else:
            return False
    
    
'''
Extending domains
'''
domain_new_boundaries = {}

for protein in plddt_list.keys(): #it's more efficient going protein by protein, instead of domain by domain
    protein_id = protein.split("-")[0]
    
    plddt_values = plddt_list[protein]
    
    protein_domains_tb = PFAM_domains[(PFAM_domains['DomainProtein'] == protein_id)]
    protein_domains_list = protein_domains_tb.set_index('index').transpose().to_dict('list') #put the datafram in dictionary (of lists) format
    
      
    for domain in protein_domains_list.keys():
        domain_id = PFAM_domains_list[domain][1]
        domain_orig_start = PFAM_domains_list[domain][2]
        domain_orig_end = PFAM_domains_list[domain][3]
        
        domain_start = PFAM_domains_list[domain][2]
        domain_end = PFAM_domains_list[domain][3]
        
        str_val = check_boundary(plddt_values,domain_start,"Start")
        while str_val == False:
            domain_start = domain_start - 1
            #print(domain_start)
            if domain_start==0: # check if the extension has reach the beggining of the protein values vector
                str_val = True
                domain_start = domain_start + 1 # if the extension already surpassed the beginnig, we go back (in) one residue 
                continue
            
            str_val = check_boundary(plddt_values,domain_start,"Start")
            
        end_val = check_boundary(plddt_values,domain_end,"End")
        end_text = "_"        
        while end_val == False:
            domain_end = domain_end + 1
            #print(domain_end)
            if domain_end>len(plddt_values): # check if the extension has reach the end of the protein values vector
                end_val = True
                domain_end = domain_end - 1 # if the window end value already surpassed the end residue, we go back (in) one residue 
                continue
            
            end_val = check_boundary(plddt_values,domain_end,"End")
        
        if PFAM_domains_list[domain][2]==domain_start:
            str_bin = 0 
        else:
            str_bin = 1 
        
        if PFAM_domains_list[domain][3]==domain_end:
            end_bin = 0  
        else:
            end_bin = 1
       
        domain_new_boundaries[domain] = [protein_id, domain_id, domain_orig_start, domain_start, domain_orig_end, domain_end, str_bin, end_bin]
        

'''
Save values in a text file
'''

print("Saving files...")


#out_file = open("HuRI_human_DMI_domains_pLDDT_extension_2.txt", "w") #name of the table
out_file = open("HuRI_human_DMI_domains_repeated_pLDDT_extension_2.txt", "w") #name of the table
out_file.write('DomainProtein\tDomainID1\tdomain_start\tdomain_pLDDT_start\tdomain_end\tdomain_pLDDT_end\tstart_ext\tend_ext\n' ) #header

# out_file = open("HumanReferenceProteome_high_pLDDTregions_70_extension.tsv", "w") #name of the table
# out_file.write('RegionProtein\tRegionNumber\tregion_start\tregion_end\tstart_ext\tend_ext\n' ) #header

for key in domain_new_boundaries.keys():
    out_file.write(f'{domain_new_boundaries[key][0]}\t{domain_new_boundaries[key][1]}\t{domain_new_boundaries[key][2]}\t{domain_new_boundaries[key][3]}\t{domain_new_boundaries[key][4]}\t{domain_new_boundaries[key][5]}\t{domain_new_boundaries[key][6]}\t{domain_new_boundaries[key][7]}\n' )
out_file.close()


'''
Check how to combine fragment valuess

Check if indexes match 

'''







