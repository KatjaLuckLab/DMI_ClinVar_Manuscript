#!/usr/bin/env Rscript

#### R script to generate titration plots for BRET experiments
# Developed by Jesús Alvarado Valverde, Luck Group
# Date created: 26/02/2026
# Date last updated: 26/02/2026
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R version 4.4.2 (2024-10-31) -- "Pile of Leaves"

# Libraries ####
library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)

# Directories ####
wd <- getwd()
myDate <- format(Sys.Date(), "%Y%m%d")
input.dir <- "Inputs"
output.dir <- "Outputs"
results.dir <- paste0("Sample_Results_",myDate,"/")

# Results dated folder
dated_Results <- dir.exists(file.path(wd, output.dir, results.dir))
if(dated_Results==FALSE){
  dir.create(file.path(wd, output.dir, results.dir))
}


# Input files - Titration data ####
# File names
fn.BRET_results <- "titration_cat_scoring_v6_with_include.tsv" # BRET result statistical analysis
fn.mutant_descriptions <- "mutant_description.csv"# Mutant information
fn.BRET_plasmid_layout <- "plate_layout.csv"  # Plate layout
fn.BRET_selection <- "titrations_filtering.csv"  # Selected BRET experiments
fn.BRET_values <- "titration_values.csv"  # BRET values
fn.BRET_fits <- "titration_fit.csv"  # BRET curve fits

# Load files
df_BRET_results <- read_tsv(file.path(wd, input.dir, fn.BRET_results), col_names = TRUE, na = c("NA","")) 
df_mutant_descriptions <- read_csv(file.path(wd, input.dir, fn.mutant_descriptions), col_names = TRUE, na = c("NA","")) 
df_BRET_plasmid_layout <- read_csv(file.path(wd, input.dir, fn.BRET_plasmid_layout), col_names = TRUE, na = c("NA","")) 
df_BRET_selection <- read_csv(file.path(wd, input.dir, fn.BRET_selection), col_names = TRUE, na = c("NA","")) 
df_BRET_values <- read_csv(file.path(wd, input.dir, fn.BRET_values), col_names = TRUE, na = c("NA","")) 
df_BRET_fits <- read_csv(file.path(wd, input.dir, fn.BRET_fits), col_names = TRUE, na = c("NA","")) 

rm(fn.BRET_results, fn.mutant_descriptions, fn.BRET_selection, fn.BRET_plasmid_layout, fn.BRET_values, fn.BRET_fits)

# Edit data frames  ####
df_BRET_selection <- df_BRET_selection %>% 
  filter(include==1 & use_fit==1) %>% 
  select(-include, -use_fit, -comment) 
df_BRET_selection$repl_id <- as.character(df_BRET_selection$repl_id)


df_BRET_plasmid_info <- df_BRET_plasmid_layout %>% 
  select(project_id, repl_id, NL_plasmid, mCit_plasmid, NL_plasmid_id, mCit_plasmid_id) %>% unique() %>%
  filter(NL_plasmid!="empty" & mCit_plasmid!="empty")
df_BRET_plasmid_info <- df_BRET_plasmid_info %>% filter(project_id %in% df_BRET_selection$project_id)

# Calculate BRET ratios
df_BRET_values$ratio <- df_BRET_values$fluo / df_BRET_values$lumi
df_BRET_values$repl_id <- as.character(df_BRET_values$repl_id)

# Add BRET values to plasmid info
df_BRET_values_info <- left_join(df_BRET_values, df_BRET_plasmid_info)

# Indicate construct names
df_BRET_values_info <- df_BRET_values_info %>% 
  mutate(NL_construct = gsub("pcDNA3.1 cmyc-NL-", "", NL_plasmid),
         mCit_construct = gsub("pcDNA3.1 mCit-His3C-", "", mCit_plasmid))


df_BRET_values_info <- df_BRET_values_info %>% 
  separate(NL_construct, c("NL_protein","NL_mutation"), sep="_", remove=F) %>% 
  mutate(NL_mutation=case_when(is.na(NL_mutation) ~ NA, !is.na(NL_mutation) ~ gsub("^.*?_", "", NL_construct) )) %>% 
  separate(mCit_construct, c("mCit_protein","mCit_mutation"), sep="_", remove=F) %>% 
  mutate(mCit_mutation=case_when(is.na(mCit_mutation) ~ NA, !is.na(mCit_mutation) ~ gsub("^.*?_", "", mCit_construct) )) %>% 
  mutate(test=case_when(is.na(NL_mutation) & is.na(mCit_mutation) ~ "WT", TRUE ~ "mutant"))


df_BRET_selection_values <- left_join(df_BRET_selection, df_BRET_values_info)
df_BRET_mutant_pairs <- df_BRET_selection_values %>% filter(test=="mutant") %>% 
  select(project_id, NL_plasmid_id, mCit_plasmid_id) %>% unique()


df_BRET_mutant_pairs <- df_BRET_results %>% 
  mutate(BRET_category=case_when(num_use_fit_WT==3 & num_use_fit_MUT==3 ~ "BRET50", 
                                     num_repl_include_MUT==3 & num_repl_include_MUT==3 ~ "maxBRET",
                                     TRUE ~ "other")) %>%
  select(project_id, NL_plasmid_id, mCit_plasmid_id, BRET_category) %>% unique()
df_BRET_mutant_pairs <- left_join(df_BRET_mutant_pairs, df_BRET_values_info %>% select(project_id, NL_plasmid_id, mCit_plasmid_id,NL_construct, mCit_construct, NL_protein, mCit_protein, test)) %>% unique()

df_BRET_mutant_pairs <- df_BRET_mutant_pairs %>% 
  mutate(gene_symbol=case_when(NL_construct==NL_protein ~ mCit_protein, mCit_construct==mCit_protein ~ NL_protein)) %>%
  mutate(mutation=case_when(NL_construct==NL_protein ~ mCit_construct, mCit_construct==mCit_protein ~ NL_construct)) %>%
  mutate(mutation=gsub("-", "_",gsub("^.*?_", "",gsub("_([0-9]{1,3}del)", "-\\1",mutation ))))

df_BRET_mutant_pairs <- left_join(df_BRET_mutant_pairs, df_mutant_descriptions  %>% mutate(canonical_mutation=gsub("-","_",canonical_mutation)))


df_BRET_mutant_pairs <- df_BRET_mutant_pairs %>% 
  mutate(mCit_construct=case_when(NL_construct==NL_protein ~ paste(gene_symbol,canonical_mutation,sep="_"), TRUE ~ mCit_construct),
         NL_construct=case_when(mCit_construct==mCit_protein ~ paste(gene_symbol,canonical_mutation,sep="_"), TRUE ~ NL_construct))


# Graphs ####
setwd(output.dir)
setwd(results.dir)
# Plot parameters
Tplot_text_size <- 18
Tplot_title_size <- 13
Tplot_width <- 3.21 # inches
Tplot_height <- 2.36 # inches

for(i in 1:nrow(df_BRET_mutant_pairs)){
  pair_info <- df_BRET_mutant_pairs[i,] %>% unlist() %>% as.vector()
  
  project <- pair_info[1]
  NL <- pair_info[2]
  mCIT <- pair_info[3]
  
  plot_data <- df_BRET_values_info %>% filter(project_id==project & NL_plasmid_id==NL & mCit_plasmid_id==mCIT & test=="mutant")
  WT_data <-  df_BRET_values_info %>% filter(project_id==project & NL_protein==pair_info[7] & mCit_protein==pair_info[8] & test=="WT")
  plot_data <- rbind(plot_data, WT_data)
  
  plot_title <- paste(gsub(" ([0-9]{1,3}del)", "-\\1", gsub("_", " ", pair_info[5])), 
                      gsub(" ([0-9]{1,3}del)", "-\\1", gsub("_", " ", pair_info[6])), sep=" - ")
  
  if(pair_info[4]=="BRET50"){
    test_label <- "BRET50"
    
    plot_fit <- df_BRET_fits %>% 
      filter(project_id==project & NL_plasmid_id==NL & mCit_plasmid_id==mCIT )%>% 
      mutate(test="mutant")
    
    WT_fit <-  df_BRET_fits %>% 
      filter(project_id==project & NL_plasmid_id==unlist(WT_data$NL_plasmid_id[1]) &
               mCit_plasmid_id==unlist(WT_data$mCit_plasmid_id[1])) %>% 
      mutate(test="WT")
    
    plot_fit <- rbind(plot_fit, WT_fit)
    
    plot_fit_data <- tibble(repl_id=character(), cBRET=double(), ratio=double(), test=character()) 
    cBRET = bret_max * expr_ratio / (bret50 + expr_ratio)
    for(i in 1:nrow(plot_fit)){
      line_info <- plot_fit[i,] %>% unlist() %>% as.vector()
      
      line_test <- line_info[10]
      replicate_id <- line_info[4]
      bret50 <- as.numeric(line_info[5])
      bretmax <-  as.numeric(line_info[7])
      
      expr_ratio_max <- plot_data %>% filter(repl_id==replicate_id & test==line_test) %>% select(ratio) %>% max()
      
      expr_ratio_line <- seq(0, expr_ratio_max, by = expr_ratio_max/100)
      cBRET_line <- bretmax * expr_ratio_line / (bret50 + expr_ratio_line)
      
      plot_fit_line <- tibble(repl_id=replicate_id, cBRET=cBRET_line, ratio=expr_ratio_line, test=line_test) 
      plot_fit_data <- rbind(plot_fit_data, plot_fit_line)
    }
    
    p <- ggplot() + geom_point(data=plot_data, size=2.2, aes(x=ratio, y = cBRET, color=test) ) + 
      geom_line(data=plot_fit_data, size=0.8, aes(x=ratio, y = cBRET, group=interaction(test, repl_id), color=test))+
      ggtitle(plot_title) +  ylab("BRET") + xlab("Protein Acc/Don") +
      annotate(geom="text", label=test_label, x=max(plot_data$ratio)*0.9, y=max(plot_data$cBRET)/3, size=3, color="red")+
      theme_classic() + labs(fill = "Pathogenicity") + scale_color_manual(values=c("#d52221","black")) +
      theme(text = element_text(size = Tplot_text_size), 
            axis.text = element_text(color = "black"), legend.position="none", 
            plot.title= element_text(size = Tplot_title_size), 
            axis.ticks = element_line(linewidth = 1, color = "black")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) + 
      scale_y_continuous(labels = function(x) number(x, accuracy = 0.1), 
                         breaks = scales::pretty_breaks(n = 3), 
                         limits=c(0,max((plot_data$cBRET))*1.1))
    
    file_name <- paste(pair_info[5], "_", pair_info[6], ".pdf", sep="")
    pdf(file_name, width = Tplot_width, height = Tplot_height)
    print(p)
    dev.off()
    
    
  } 
  else if(pair_info[4]=="maxBRET"){
    test_label <- "maxBRET"
    
    p <- ggplot(plot_data, aes(x=ratio, y = cBRET, color=test) ) + 
      geom_point(size=2.2) +
      ggtitle(plot_title) +  ylab("BRET") + xlab("Protein Acc/Don") +
      annotate(geom="text", label=test_label, x=max(plot_data$ratio)*0.9, y=max(plot_data$cBRET)/3, size=3, color="red")+
      theme_classic() + labs(fill = "Pathogenicity") + scale_color_manual(values=c("#d52221","black")) +
      theme(text = element_text(size = Tplot_text_size), 
            axis.text = element_text(color = "black"), 
            legend.position="none", 
            plot.title= element_text(size = Tplot_title_size), 
            axis.ticks = element_line(linewidth = 1, color = "black")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 2)) + 
      scale_y_continuous(labels = function(x) number(x, accuracy = 0.1), 
                         breaks = scales::pretty_breaks(n = 3), 
                         limits=c(0, max((plot_data$cBRET))*1.1)
                         ) 
    
    file_name <- paste(pair_info[5], "_",pair_info[6], ".pdf", sep="")
    pdf(file_name, width = Tplot_width, height = Tplot_height)
    print(p)
    dev.off()
    
    
  } 
}






