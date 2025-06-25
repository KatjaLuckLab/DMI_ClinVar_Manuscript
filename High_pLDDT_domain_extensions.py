#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:25:09 2024

@author: jesusav
"""

import os
import re
import pandas as pd


'''
Make a list of all proteins in the proteome

'''

os.chdir('/Volumes/imb-luckgr/imb-luckgr2/projects/AlphaFold/AFDB_HumanReferenceProteome/')
#os.chdir('.')

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
Region parameters
'''
region_minL = 20 # Residue window in which the pLDDT values would be averaged, checked and compared
pLDDT_region_thrs = 70 # pLDDT value that the domain window has to overcome in order to be extended (to be kept)



'''
Find regions of high pLDDT values and store them in a protein region dictionary
'''

plddt_regions = {}

for protein in plddt_list.keys():
    residue_n = 1
    residue_window_init = 0  # start/initialise a high pLDDT value window
    region_start = None
    region_end = None
    region_boundaries = {}
    
    for residue in plddt_list[protein]: # loop throough all pLDDT values    
        if residue >= pLDDT_region_thrs: # check if the value of a given residue if above the threshold
            residue_window_init += 1 # start/initialise a high pLDDT value window 
            
            if residue_window_init == 2: #if there are two consecutive high values we define a region start residue
                region_start = residue_n - 1 # the region would have started on the previous residue so we rest one
                # we do this here so that we don't have to update it later on, as the start would be independent
                
            if residue_window_init >= region_minL: # if the window reaches a length of 20 then we define a region end residue
                region_end = residue_n 
                region_boundaries[region_start] = region_end # we then create/update a dictionary of region of the current protein
                #with the region start as a key and the end as the value to be updated if the window keeps on extending 
                
        else: #if the residue pLDDT value is lower than the threshold 
            residue_window_init = 0 # we restart the window
            region_end = None # we restart the window end
            
        residue_n += 1 # increase the residue counter
    
    # after looping through the pLDDT values
    if len(region_boundaries) > 0:  # we loop through the dictionary of regions per protein, if there is one with region values
        region_n = 1 # make a counter for the regions in the dictionary
        
        for region in region_boundaries.keys(): 
            region_key = protein + "_R" + str(region_n) # create a region key comprised by the protein Id and teh region number
            plddt_regions[region_key] = [region, region_boundaries[region]] # store start and end of the region in the general regions dictionary
            region_n += 1

del region, region_n, region_key
del residue, residue_n, residue_window_init, region_start, region_end,
del protein, region_boundaries



'''
Save region boundaries in a text file
'''

print("Saving files...")


out_file = open("HumanReferenceProteome_high_pLDDTregions_70.tsv", "w") #name of the table
out_file.write('Region_protein\tRegion_code\tregion_start\tregion_end\n' ) #header
for region_key in plddt_regions.keys():
    protein = region_key.split("_")
    out_file.write(f'{protein[0]}\t{protein[1]}\t{plddt_regions[region_key][0]}\t{plddt_regions[region_key][1]}\n' )
out_file.close()

del region_key, protein








