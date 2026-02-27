##########################################################
# File: BRET_heatmaps_DMI_functions.R
# Author: Jesus_Alvarado_Valverde
# Date: 2026-02-27
# Description: Functions to create networks tables and graphs 
#               using the data format of DMIs and DMI AF models
##########################################################


# ==========================
# Dependencies
# ==========================

library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(colorspace)
library(ggh4x)

# ==========================
# General variables
# ==========================
# Plot palettes
colors_BRET <- c("grey40", "green", "magenta")
colors_Var <-c("#56B4E9","#E69F00","#d52221","#999999")

# Heatmap parameters
HM_parameters <- c(Hplot_text_size=20,  
                   legend_t_size=15, 
                   legend_t_vjust=3, 
                   legend_tx_size <- 13)
heatmap_pltt='RdBu'

heatmap_map <- c(
  tstat_bret50 = "Effect size\nBRET50",
  tstat_cBRET = "Effect size\nmaxBRET",
  BRET50_Log2_FChng = "Log2 FC\nBRET50",
  maxBRET_Log2_FChng = "Log2 FC\nmaxBRET"
)

# ==========================
# Create heat map based on BRET analysis of domain mutations among DMI PPIs
# ==========================

#' This function takes a DMI BRET results table and save a PDF file with a heatmap for the domain mutations
#'
#' @param domain_BRET A tibble with data on the BRET test of domain mutations
#' @param proteinSymbol A character variable indicating the name of the protein group
#' @param heatmap_type A character variable indicating the name of the type of the heat map between effect sizes or log2FC, between BRET50 and maxBRET
#' @return No returns, the function will automatically generate a PDF file in the current working directoy
#' @examples
#' BRETheatmap_domain(domain_BRET, "CTBP1", "maxBRET_Log2_FChng") # Returns ...
#' 
BRETheatmap_domain <- function(domain_BRET, proteinSymbol, heatmap_type){
  # Define value limits and significance labels based on the type of domain heatmap
  if(heatmap_type=="tstat_bret50"){
    heatmap_limit <- 20
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="tstat_cBRET"){
    heatmap_limit <- 20
    sig_label <- "sig_maxBRET"
  }else if (heatmap_type=="BRET50_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="maxBRET_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_maxBRET"
  }
  
  df_heatmap <- domain_BRET %>% filter(NL_symbol == proteinSymbol)
  heatmap_plot <- ggplot(df_heatmap, aes(fct_reorder(canonical_mutation, residue), mCit_symbol, fill = !!as.symbol(heatmap_type))) + 
    geom_tile(color = "white", lwd = 1.2) + 
    geom_text(aes(label = !!as.symbol(sig_label)), size = 7) +
    scale_fill_continuous_divergingx(palette = heatmap_pltt, rev = TRUE, mid = 0, 
                                     limits = c(-heatmap_limit, heatmap_limit), na.value = "grey80")+
    xlab(paste(proteinSymbol, "mutants")) + theme_minimal() + 
    labs(fill = heatmap_map[heatmap_type]) +
    theme(text = element_text(size =  HM_parameters["Hplot_text_size"]), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_text(color = "black"), 
          axis.ticks = element_line(linewidth = 0), 
          axis.title.y = element_blank(), 
          axis.title.x = element_text(hjust = 0.3),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), 
          plot.title = element_blank(),
          legend.key.height = unit(0.6, 'cm'), 
          legend.key.width = unit(0.25, 'cm'), 
          legend.title = element_text(size =  HM_parameters["legend_t_size"], vjust = HM_parameters["legend_t_vjust"]), 
          legend.text = element_text(size = HM_parameters["legend_tx_size"])) +
    scale_y_discrete(limits=rev)  + 
    facet_grid2(category~., scales = "free_y", switch = "y", space = "free_y",
                strip = strip_themed(text_y = element_text(size=0),
                                     background_y = elem_list_rect( fill =colors_BRET )))
  # SAVE PDF
  Hplot_width <- 8 + nrow(df_heatmap %>% select(canonical_mutation) %>% unique()) # 8 minimum width + the number of mutations
  Hplot_height <- 6.5 + 0.5*nrow(df_heatmap %>% select(mCit_symbol) %>% unique()) # 6 minimum height + the number of proteins
  Hplot_filename <- paste("DomainValidation_", proteinSymbol, "_", heatmap_type, ".pdf", sep="")
  ggsave(heatmap_plot, filename = Hplot_filename, device = "pdf", height = Hplot_height, width = Hplot_width, units = "cm")
  
}

# ==========================
# Create heat map based on BRET analysis of motif mutations among DMI PPIs
# ==========================

#' This function takes a DMI BRET results table and save a PDF file with a heatmap for the motif deletions
#'
#' @param motif_BRET A tibble with data on the BRET test of motif mutations
#' @param proteinSymbol A character variable indicating the name of the protein group
#' @param heatmap_type A character variable indicating the name of the type of the heat map between effect sizes or log2FC, between BRET50 and maxBRET
#' @return No returns, the function will automatically generate a PDF file in the current working directoy
#' @examples
#' BRETheatmap_motif(motif_BRET, "CTBP1", "maxBRET_Log2_FChng") # Returns ...
#' 
BRETheatmap_motif <- function(motif_BRET, proteinSymbol, heatmap_type){
  # Define value limits and significance labels based on the type of domain heatmap
  if(heatmap_type=="tstat_bret50"){
    heatmap_limit <- 20
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="tstat_cBRET"){
    heatmap_limit <- 20
    sig_label <- "sig_maxBRET"
  }else if (heatmap_type=="BRET50_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="maxBRET_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_maxBRET"
  }
  
  df_heatmap <- motif_BRET %>% filter(NL_symbol==proteinSymbol) %>% mutate(label=paste(gene_symbol,canonical_mutation, sep="\n"))
  heatmap_plot <- ggplot(df_heatmap, aes(NL_symbol, fct_reorder2(label, as.numeric(residue), desc(mCit_symbol)), fill=!!as.symbol(heatmap_type))) + 
    geom_tile(color = "white", lwd = 1.2) + 
    geom_text(aes(label=!!as.symbol(sig_label)), size=7)  +
    scale_fill_continuous_divergingx(palette = heatmap_pltt, rev = TRUE, mid = 0, 
                                     limits = c(-heatmap_limit, heatmap_limit), na.value = "grey80")+
    ylab("Motif deletions") + theme_minimal() + 
    labs(fill = heatmap_map[heatmap_type]) +
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(), 
          axis.text = element_text(color = "black"), 
          panel.grid.minor = element_blank(), 
          axis.ticks = element_line(size = 0), 
          axis.title.x=element_blank(), 
          plot.title=element_text(size=12, hjust=0.95),
          legend.key.height= unit(0.5, 'cm'), 
          legend.key.width= unit(0.25, 'cm'), 
          legend.title = element_text(size = HM_parameters["legend_t_size"], vjust = HM_parameters["legend_t_vjust"]), 
          legend.text = element_text(size = HM_parameters["legend_tx_size"])) +
    scale_y_discrete(limits=rev)+ 
    facet_grid2(category~., scales = "free_y", switch = "y",space = "free_y",
                strip = strip_themed(text_y = element_text(size=0),
                                     background_y = elem_list_rect( fill =colors_BRET[2:3] ))) 
  
  # SAVE PDF
  Hplot_width <- 8 + nrow(df_heatmap %>% select(mCit_symbol) %>% unique()) # 8 minimum width + the number of mutations
  Hplot_height <- 4.5 + nrow(df_heatmap %>% select(canonical_mutation) %>% unique()) # 6 minimum height + the number of proteins
  Hplot_filename <- paste("MotifValidation_", proteinSymbol, "_", heatmap_type, ".pdf", sep="")
  ggsave(heatmap_plot, filename = Hplot_filename, device = "pdf", height = Hplot_height, width = Hplot_width, units = "cm")
  
}

# ==========================
# Create heat map based on BRET analysis of domain variants among DMI PPIs
# ==========================

#' This function takes a DMI BRET results table and save a PDF file with a heatmap for the domain variants
#'
#' @param variants_BRET A tibble with data on the BRET test of protein variants
#' @param proteinSymbol A character variable indicating the name of the protein group
#' @param heatmap_type A character variable indicating the name of the type of the heat map between effect sizes or log2FC, between BRET50 and maxBRET
#' @return No returns, the function will automatically generate a PDF file in the current working directoy
#' @examples
#' BRETheatmap_VarDomain(variants_BRET, "CTBP1", "maxBRET_Log2_FChng") # Returns ...
#' 
BRETheatmap_VarDomain <- function(variants_BRET, proteinSymbol, heatmap_type){
  # Define value limits and significance labels based on the type of domain heatmap
  if(heatmap_type=="tstat_bret50"){
    heatmap_limit <- 20
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="tstat_cBRET"){
    heatmap_limit <- 20
    sig_label <- "sig_maxBRET"
  }else if (heatmap_type=="BRET50_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="maxBRET_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_maxBRET"
  }
  
  df_heatmap <- variants_BRET %>% filter(NL_symbol == proteinSymbol & domain_variant==TRUE)
  
  variant_path <- df_heatmap %>% select(ClinVar) %>% unique() %>% unlist()
  variant_colors <- colors_Var[3:4] 
  if(length(variant_path)<2){
    if(variant_path=="Pathogenic"){
      variant_colors <- colors_Var[3]
    } else if (variant_path=="Uncertain"){
      variant_colors <- colors_Var[4]
    }
  } 
  
  heatmap_plot <- ggplot(df_heatmap, aes(fct_reorder(canonical_mutation, residue), mCit_symbol, fill=!!as.symbol(heatmap_type))) + 
    geom_tile(color = "white", lwd = 1.2) + 
    geom_text(aes(label = !!as.symbol(sig_label)), size = 7) +
    scale_fill_continuous_divergingx(palette = heatmap_pltt, rev = TRUE, mid = 0, 
                                     limits = c(-heatmap_limit, heatmap_limit), na.value = "grey80")+
    xlab(paste(proteinSymbol, "variants")) + theme_minimal() + 
    labs(fill = heatmap_map[heatmap_type]) +
    theme(text = element_text(size =  HM_parameters["Hplot_text_size"]), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text = element_text(color = "black"), 
          axis.ticks = element_line(linewidth = 0), 
          axis.title.y = element_blank(), 
          axis.title.x = element_text(hjust = 0.3),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), 
          plot.title = element_blank(),
          legend.key.height = unit(0.6, 'cm'), 
          legend.key.width = unit(0.25, 'cm'), 
          legend.title = element_text(size =  HM_parameters["legend_t_size"], vjust = HM_parameters["legend_t_vjust"]), 
          legend.text = element_text(size = HM_parameters["legend_tx_size"])) +
    scale_y_discrete(limits=rev)  + 
    facet_grid2(.~ClinVar, scales = "free_x", switch = "x", space = "free_x",
                strip = strip_themed(text_x = element_text(size=0),
                                     background_x = elem_list_rect( fill = variant_colors )))
  # SAVE PDF
  Hplot_width <- 8 + nrow(df_heatmap %>% select(canonical_mutation) %>% unique()) # 8 minimum width + the number of mutations
  Hplot_height <- 6.5 + 0.5*nrow(df_heatmap %>% select(mCit_symbol) %>% unique()) # 6 minimum height + the number of proteins
  Hplot_filename <- paste("Variant_Domain_", proteinSymbol, "_", heatmap_type, ".pdf", sep="")
  ggsave(heatmap_plot, filename = Hplot_filename, device = "pdf", height = Hplot_height, width = Hplot_width, units = "cm")
  
}


# ==========================
# Create heat map based on BRET analysis of motif variants among DMI PPIs
# ==========================

#' This function takes a DMI BRET results table and save a PDF file with a heatmap for the motif variants
#'
#' @param variants_BRET A tibble with data on the BRET test of protein variants
#' @param proteinSymbol A character variable indicating the name of the protein group
#' @param heatmap_type A character variable indicating the name of the type of the heat map between effect sizes or log2FC, between BRET50 and maxBRET
#' @return No returns, the function will automatically generate a PDF file in the current working directoy
#' @examples
#' BRETheatmap_VarMotif(variants_BRET, "CTBP1", "maxBRET_Log2_FChng") # Returns ...
#' 
BRETheatmap_VarMotif <- function(variants_BRET, proteinSymbol, heatmap_type){
  # Define value limits and significance labels based on the type of domain heatmap
  if(heatmap_type=="tstat_bret50"){
    heatmap_limit <- 20
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="tstat_cBRET"){
    heatmap_limit <- 20
    sig_label <- "sig_maxBRET"
  }else if (heatmap_type=="BRET50_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_BRET50"
  }else if (heatmap_type=="maxBRET_Log2_FChng"){
    heatmap_limit <- 4
    sig_label <- "sig_maxBRET"
  }
  
  df_heatmap <- variants_BRET %>% filter(NL_symbol == proteinSymbol & motif_variant==TRUE) %>% mutate(label=paste(gene_symbol,canonical_mutation, sep="\n"))
  
  variant_path <- df_heatmap %>% select(ClinVar) %>% unique() %>% unlist()
  
  variant_colors <- colors_Var[3:4] 
  if(length(variant_path)<2){
    if(variant_path=="Pathogenic"){
      variant_colors <- colors_Var[3]
    } else if (variant_path=="Uncertain"){
      variant_colors <- colors_Var[4]
    }
  } 
  
  heatmap_plot <- ggplot(df_heatmap, aes(NL_symbol, fct_reorder2(label, residue, mCit_symbol), fill=!!as.symbol(heatmap_type))) + 
    geom_tile(color = "white", lwd = 1.2) + 
    geom_text(aes(label=!!as.symbol(sig_label)), size=7)  +
    scale_fill_continuous_divergingx(palette = heatmap_pltt, rev = TRUE, mid = 0, 
                                     limits = c(-heatmap_limit, heatmap_limit), na.value = "grey80")+
    ylab("Motif variants") + theme_minimal() + 
    labs(fill = heatmap_map[heatmap_type]) +
    theme(text = element_text(size = 20), 
          panel.grid.major = element_blank(), 
          axis.text = element_text(color = "black"), 
          panel.grid.minor = element_blank(), 
          axis.ticks = element_line(size = 0), 
          axis.title.x=element_blank(), 
          plot.title=element_blank(),
          legend.key.height= unit(0.5, 'cm'), 
          legend.key.width= unit(0.25, 'cm'), 
          legend.title = element_text(size = HM_parameters["legend_t_size"], vjust = HM_parameters["legend_t_vjust"]), 
          legend.text = element_text(size = HM_parameters["legend_tx_size"])) +
    scale_y_discrete(limits=rev)+ 
    facet_grid2(ClinVar~., scales = "free_y", switch = "y", space = "free_y",
                strip = strip_themed(text_y = element_text(size=0),
                                     background_y = elem_list_rect( fill=variant_colors))) 
  # SAVE PDF
  
  Hplot_width <- 8 + nrow(df_heatmap %>% select(mCit_symbol) %>% unique()) # 8 minimum width + the number of mutations
  Hplot_height <- 4.5 + nrow(df_heatmap %>% select(canonical_mutation) %>% unique()) # 6 minimum height + the number of proteins
  Hplot_filename <- paste("Variant_Motif_", proteinSymbol, "_", heatmap_type, ".pdf", sep="")
  ggsave(heatmap_plot, filename = Hplot_filename, device = "pdf", height = Hplot_height, width = Hplot_width, units = "cm")
  
}


