#!/usr/bin/env Rscript

#### R script to generate Heatmaps plots for BRET experiments
# Developed by Jesús Alvarado Valverde, Luck Group, IMB
# Date created: 26/02/2026
# Date last updated: 27/02/2026
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R version 4.4.2 (2024-10-31) -- "Pile of Leaves"

# Libraries ####
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(colorspace)
library(ggh4x)

# Directories ####
wd <- getwd()
myDate <- format(Sys.Date(), "%Y%m%d")
input.dir <- "Inputs"
output.dir <- "Outputs"
results.dir <- paste0("Heatmaps_",myDate,"/")

# Results dated folder
dated_Results <- dir.exists(file.path(wd, output.dir, results.dir))
if(dated_Results==FALSE){
  dir.create(file.path(wd, output.dir, results.dir))
}


# Input files - BRET data ####
# File names
fn.BRET_results <- "titration_cat_scoring_v6_with_include.tsv" # BRET result statistical analysis
fn.mutant_descriptions <- "mutant_description.csv"# Mutant information
fn.BRET_mutants <- "tested_disease_mutations_structure_computations_extended.txt" # Information about the variant pathogenicity

# Load files
df_BRET_results <- read_tsv(file.path(wd, input.dir, fn.BRET_results), col_names = TRUE, na = c("NA","")) 
df_mutant_descriptions <- read_csv(file.path(wd, input.dir, fn.mutant_descriptions), col_names = TRUE, na = c("NA","")) 
df_BRET_mutants <- read_tsv(file.path(wd, input.dir, fn.BRET_mutants), col_names = TRUE, na = c("NA","")) 

rm(fn.BRET_results, fn.mutant_descriptions, fn.BRET_mutants)

# Edit data frames  ####
df_BRET_results <- df_BRET_results %>% # Correct and exclude cases 
  mutate(mut_type=case_when(mCit_plasmid=="pcDNA3.1 mCit-His3C-MYD88_134_138del" ~ "Motif validation", TRUE ~ mut_type))

df_BRET_results <- df_BRET_results %>% 
  mutate(include=case_when(NL_plasmid=="pcDNA3.1 cmyc-NL-IQCB1_N273Y" &  mCit_symbol=="MNS1" & project_id=="Lu153r01" ~ "1",
                           NL_plasmid=="pcDNA3.1 cmyc-NL-SPOP_S119L" &  mCit_symbol=="MYD88" & project_id=="Lu171r01" ~ "0",
                           NL_plasmid=="pcDNA3.1 cmyc-NL-SPOP_R70T" &  mCit_symbol=="RXRB" & project_id=="Lu169r01" ~ "0" , TRUE ~ include))

df_BRET_results <- df_BRET_results %>% filter(include=="1") # Select PPI pairs to include in the plotting

df_BRET_results <- df_BRET_results %>% # Label type of experimental category
  mutate(exp_category=case_when(mut_type=="Domain validation" | mut_type=="Motif validation" ~ "interface", TRUE ~ "variant"))

#Significance labels
df_BRET_results <- df_BRET_results %>%  
  mutate(tstat_bret50 = as.double(gsub("\\,","\\.",tstat_bret50))) %>% # Convert effect size and pvalue formats from "," to "."
  mutate(tstat_cBRET = as.double(gsub("\\,","\\.",tstat_cBRET)) ) %>% 
  mutate(pval_bret50 = as.double(gsub("\\,","\\.",pval_bret50))) %>% 
  mutate(pval_cBRET = as.double(gsub("\\,","\\.",pval_cBRET)) ) %>% 
  mutate(sig_BRET50=case_when(pval_bret50 > 0.05 ~ "",  # Assign significance level labels
                              pval_bret50 <= 0.05 & pval_bret50 > 0.01 ~ "*", 
                              pval_bret50 <= 0.01 & pval_bret50 > 0.001 ~ "**",
                              pval_bret50 <= 0.001 ~ "***")) %>%
  mutate(sig_maxBRET=case_when(pval_cBRET>0.05 ~ "",  
                               pval_cBRET <= 0.05 & pval_cBRET > 0.01 ~ "*", 
                               pval_cBRET <= 0.01 & pval_cBRET > 0.001 ~ "**",
                               pval_cBRET <= 0.001 ~ "***"))

df_BRET_results <- df_BRET_results %>% # Label cases when fit for BRET50 was successful
  mutate(BRET_category=case_when(num_use_fit_WT == 3 & num_use_fit_MUT == 3 ~ "BRET50", 
                                 num_repl_include_MUT == 3 & num_repl_include_MUT == 3 ~ "maxBRET",
                                 TRUE ~ "other"))

### Fold change calculations
df_BRET_results <- df_BRET_results  %>%  
  mutate(avg_bret50_MUT = as.double(gsub("\\,","\\.",avg_bret50_MUT)) ) %>% 
  mutate(avg_bret50_WT = as.double(gsub("\\,","\\.",avg_bret50_WT)) )%>% 
  mutate(BRET50_Log2_FChng = log2(avg_bret50_MUT/avg_bret50_WT)) 

df_BRET_results <- df_BRET_results  %>%  
  mutate(avg_cBRET_WT = as.double(gsub("\\,","\\.",avg_cBRET_WT)) ) %>% 
  mutate(avg_cBRET_MUT = as.double(gsub("\\,","\\.",avg_cBRET_MUT)) )%>% 
  mutate(maxBRET_Log2_FChng = log2(avg_cBRET_MUT/avg_cBRET_WT)) 

### Mutant descriptions
df_mutant_descriptions <- df_mutant_descriptions %>% 
  select(-template_plasmid_name, -description, -success, -alpha_missense_pred) #Remove unused variables

### Domain validation data
domain_BRET_categories <- df_BRET_results %>% select(mCit_symbol) %>% unique()
domain_BRET_categories$category <- as.factor(c("candidate", "positive", "negative", 
                                               "positive", "candidate", "negative",
                                               "candidate","negative",
                                               "positive", "candidate",
                                               "candidate", "candidate","positive","negative","negative"))


domain_BRET <- df_BRET_results %>% 
  filter(mut_type=="Domain validation") %>% 
  mutate(mutant=gsub("^pc.*?NL-", "", NL_plasmid)) %>% # Mutant labels
  mutate(residue=as.numeric(gsub("[A-Z]", "", gsub("^.*?_", "", mutant)))) %>% # Mutant index
  mutate(tstat_bret50=case_when(tstat_bret50 < -20 ~ -20, tstat_bret50 > 20 ~ 20, TRUE ~ tstat_bret50), 
         tstat_cBRET=case_when(tstat_cBRET < -20 ~ -20, tstat_cBRET > 20 ~ 20, TRUE ~ tstat_cBRET)) # Capping values

domain_BRET <- domain_BRET %>% 
  select(mCit_plasmid,NL_symbol, mCit_symbol, pval_bret50, tstat_bret50, pval_cBRET, tstat_cBRET, sig_BRET50, sig_maxBRET, BRET_category, mutant, residue, BRET50_Log2_FChng, maxBRET_Log2_FChng)

domain_BRET <- left_join(domain_BRET, domain_BRET_categories)
domain_BRET$category <- factor(domain_BRET$category, levels=c("negative", "positive", "candidate"))
domain_BRET <- domain_BRET %>% 
  mutate(mutant=case_when(mutant=="IQCB1" ~ gsub("pcDNA3.1 mCit-His3C-","",mCit_plasmid), TRUE~ mutant))%>% 
  separate(mutant, c("gene_symbol","mutation"), sep="_", remove=F)
domain_BRET <- left_join(domain_BRET, df_mutant_descriptions)


### Motif validation data
motif_BRET <- df_BRET_results %>% 
  filter(mut_type=="Motif validation") %>% 
  mutate(mutant=gsub("^pc.*?His3C-", "", mCit_plasmid)) %>%  # Mutant labels
  mutate(residue=as.numeric(gsub("[A-Z]", "", gsub("_.*?$", "",gsub("^.*?_", "", mutant))))) %>% # Mutant index
  mutate(tstat_bret50=case_when(tstat_bret50 < -20 ~ -20, tstat_bret50 > 20 ~ 20, TRUE ~ tstat_bret50), 
         tstat_cBRET=case_when(tstat_cBRET < -20 ~ -20, tstat_cBRET > 20 ~ 20, TRUE ~ tstat_cBRET)) # Capping values

motif_BRET <- motif_BRET %>% 
  select(NL_plasmid, NL_symbol, mCit_symbol, pval_bret50, tstat_bret50, pval_cBRET, tstat_cBRET, sig_BRET50, sig_maxBRET, BRET_category, mutant, residue, BRET50_Log2_FChng, maxBRET_Log2_FChng)

motif_BRET <- left_join(motif_BRET, domain_BRET_categories)
motif_BRET$category <- factor(motif_BRET$category, levels=c("negative", "positive", "candidate"))

motif_BRET <- motif_BRET  %>% 
  mutate(mutant=case_when(NL_symbol=="IQCB1" ~ gsub("pcDNA3.1 cmyc-NL-","",NL_plasmid), TRUE~ mutant))%>% 
  mutate(mutant=gsub("(^.*?_[0-9]{1,3})_","\\1-",mutant)) %>% 
  separate(mutant, c("gene_symbol","mutation"), sep="_", remove=F)
motif_BRET <- left_join(motif_BRET, df_mutant_descriptions %>% mutate(mutation=gsub("_","-",mutation), canonical_mutation=gsub("_","-",canonical_mutation)))
motif_BRET <- motif_BRET %>% 
  mutate(canonical_mutation=case_when(mutant=="MYD88_134-138del" ~ "134-138del", TRUE ~ canonical_mutation), 
         uniprot_id=case_when(mutant=="MYD88_134-138del" ~ "Q99836", TRUE ~ uniprot_id))

### Mutants testing data
mutants_BRET <- df_BRET_results %>% 
  filter(exp_category=="variant") %>% 
  mutate(domain_variant=str_detect(NL_plasmid, "_[A-Z][0-9]{1,4}[A-Z]"), motif_variant=str_detect(mCit_plasmid, "_[A-Z][0-9]{1,4}[A-Z]")) %>% # mutant labels
  mutate(mutant=case_when(domain_variant==TRUE ~ gsub("^pc.*?NL-", "", NL_plasmid), 
                          motif_variant==TRUE ~ gsub("^pc.*?His3C-", "", mCit_plasmid) )) %>% 
  mutate(mutant=case_when(domain_variant==TRUE ~ gsub("_", " ",mutant), motif_variant==TRUE ~ gsub("_", " ",mutant))) %>% 
  mutate(residue=as.numeric(gsub("[A-Z]", "", gsub("^.*? ", "", mutant)))) %>% # Mutant index
  mutate(tstat_bret50=case_when(tstat_bret50 < -20 ~ -20, tstat_bret50 > 20 ~ 20, TRUE ~ tstat_bret50), 
         tstat_cBRET=case_when(tstat_cBRET < -20 ~ -20, tstat_cBRET > 20 ~ 20, TRUE ~ tstat_cBRET)) # Capping values

mutants_BRET <- mutants_BRET %>% 
  select(NL_plasmid, mCit_plasmid,NL_symbol, mCit_symbol, pval_bret50, tstat_bret50, pval_cBRET, tstat_cBRET, sig_BRET50, sig_maxBRET, BRET_category, motif_variant, domain_variant, mutant, residue, BRET50_Log2_FChng, maxBRET_Log2_FChng)

mutants_BRET_categories <- mutants_BRET %>% select(mutant) %>% unique()

#BRET mutants
Clinvar_mutants <- df_BRET_mutants %>% 
  mutate(mutant=case_when(mutation_location=="domain" ~ paste(domain_partner,aa_change, sep=" "), 
                          TRUE ~ paste(motif_partner,aa_change, sep=" "))) %>%
  select(mutant, ClinVar) %>% unique()
mutants_BRET_categories <- left_join(mutants_BRET_categories, Clinvar_mutants)
mutants_BRET <- left_join(mutants_BRET, mutants_BRET_categories)
mutants_BRET <- left_join(mutants_BRET, domain_BRET_categories)
mutants_BRET$category <- factor(mutants_BRET$category, levels=c("negative", "positive", "candidate"))
mutants_BRET <- mutants_BRET %>% mutate(ClinVar=case_when(mutant=="IKZF1 S41L" ~ "Uncertain", 
                                                          mutant=="SPOP Y87C" ~ "Uncertain", TRUE ~ ClinVar))
mutants_BRET <- mutants_BRET %>% separate(mutant, c("gene_symbol","mutation"), sep=" ", remove=F)
mutants_BRET <- left_join(mutants_BRET, df_mutant_descriptions)


## Ploting ####
setwd("/Users/jesusav/Desktop")
source("BRET_heatmaps_DMI_functions.R")

setwd(file.path(output.dir, results.dir)) 

heatmap_types <- c("tstat_bret50", "tstat_cBRET", "BRET50_Log2_FChng", "maxBRET_Log2_FChng")
domainSymbols <- c("CTBP1", "SPOP", "WWOX", "PPP3CA")

# Domain graphs 
for(i in 1:length(domainSymbols)){
    BRETheatmap_domain(domain_BRET, domainSymbols[i], heatmap_types[4])
}
BRETheatmap_motif(domain_BRET, "IQCB1", heatmap_types[4])

# Motif graphs
for(i in 1:length(domainSymbols)){
  BRETheatmap_motif(motif_BRET, domainSymbols[i], heatmap_types[4])
}
BRETheatmap_domain(motif_BRET, "IQCB1", heatmap_types[4])

# Variant graphs
VarDomain_groups <- mutants_BRET %>% filter(domain_variant==TRUE) %>% select(NL_symbol) %>% unique() %>% unlist()
for(i in 1:length(VarDomain_groups)){
  BRETheatmap_VarDomain(mutants_BRET, VarDomain_groups[i], heatmap_types[4])
}

VarMotif_groups <- mutants_BRET %>% filter(motif_variant==TRUE) %>% select(NL_symbol) %>% unique() %>% unlist()
for(i in 1:length(VarMotif_groups)){
  BRETheatmap_VarMotif(mutants_BRET, VarMotif_groups[i], heatmap_types[4])
}



